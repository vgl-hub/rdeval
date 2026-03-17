#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cstdlib>

#include "log.h"
#include "global.h"
#include "bed.h"
#include "struct.h"
#include "gfa-lines.h"
#include "uid-generator.h"
#include "gfa.h"
#include "functions.h" // global functions
#include "stream-obj.h"

#include "cifi.h"
#include "bam-utils.h"

// Find enzyme info from name
EnzymeInfo get_enzyme(const std::string& name) {
	// Adjust cut_offset if you want exact biological cleavage positions.
	if (name == "HindIII")
		return EnzymeInfo{"HindIII", "AAGCTT", 1};
	else if (name == "NlaIII")
		return EnzymeInfo{"NlaIII", "CATG", 4};
	else if (name == "DpnII")
		return EnzymeInfo{"DpnII", "GATC", 0};
	else
		throw std::runtime_error("Unsupported restriction enzyme: " + name);
}

std::vector<int> find_cut_positions(const InRead& read, const EnzymeInfo& enz) {
	if (!read.inSequence) {
		throw std::runtime_error("InRead::inSequence is null in find_cut_positions");
	}

	const std::string& seq = *read.inSequence;
	std::vector<int> cuts;

	const std::string& motif = enz.motif;
	const int m = static_cast<int>(motif.size());
	const int n = static_cast<int>(seq.size());

	for (int i = 0; i <= n - m; ++i) {
		bool match = true;
		for (int j = 0; j < m; ++j) {
			if (seq[i + j] != motif[j]) {
				match = false;
				break;
			}
		}
		if (match) {
			int cut_pos = i + enz.cut_offset;
			if (cut_pos > 0 && cut_pos < n) {
				cuts.push_back(cut_pos);
			}
		}
	}
	return cuts;
}

void digest_read(const InRead& read, const EnzymeInfo& enz, std::vector<std::string>& frag_seqs, std::vector<std::string>& frag_quals, std::vector<int>* frag_lengths) {
	if (!read.inSequence) {
		throw std::runtime_error("InRead::inSequence is null in digest_read");
	}
	if (!read.inSequenceQuality) {
		throw std::runtime_error("InRead::inSequenceQuality is null in digest_read");
	}

	const std::string& seq  = *read.inSequence;
	const std::string& qual = *read.inSequenceQuality;

	frag_seqs.clear();
	frag_quals.clear();
	if (frag_lengths) frag_lengths->clear();

	const int n = static_cast<int>(seq.size());
	if (n == 0) return;

	if (qual.size() != seq.size()) {
		throw std::runtime_error("Sequence and quality length mismatch in digest_read");
	}

	std::vector<int> cut_positions = find_cut_positions(read, enz);
	std::vector<int> boundaries;
	boundaries.reserve(cut_positions.size() + 2);

	boundaries.push_back(0);
	for (int c : cut_positions) {
		if (c > 0 && c < n)
			boundaries.push_back(c);
	}
	boundaries.push_back(n);

	for (std::size_t i = 0; i + 1 < boundaries.size(); ++i) {
		int start = boundaries[i];
		int end   = boundaries[i + 1];
		if (end > start) {
			frag_seqs.push_back(seq.substr(start, end - start));
			frag_quals.push_back(qual.substr(start, end - start));
			if (frag_lengths)
				frag_lengths->push_back(end - start);
		}
	}
}

void chop_read_into_bambatch(const InRead& read, const EnzymeInfo& enz, BamBatch& batch) {
	if (!read.inSequenceQuality) {
		throw std::runtime_error("InRead::inSequenceQuality is null in chop_read_into_bambatch");
	}

	std::vector<std::string> frag_seqs;
	std::vector<std::string> frag_quals;
	digest_read(read, enz, frag_seqs, frag_quals, nullptr);

	// Ensure enough capacity in one go
	bambatch_ensure_capacity(batch, frag_seqs.size());

	for (std::size_t i = 0; i < frag_seqs.size(); ++i) {
		std::string frag_id = read.seqHeader + "_frag" + std::to_string(i);
		bambatch_append_unmapped(batch, frag_id, frag_seqs[i], frag_quals[i], 0);
	}
}

void build_pe_bambatches(const InRead& read, const EnzymeInfo& enz,BamBatch& R1_batch, BamBatch& R2_batch) {
	if (!read.inSequenceQuality) {
		throw std::runtime_error("InRead::inSequenceQuality is null in build_pe_bambatches");
	}

	std::vector<std::string> frag_seqs;
	std::vector<std::string> frag_quals;
	digest_read(read, enz, frag_seqs, frag_quals, nullptr);

	const std::size_t n = frag_seqs.size();
	if (n == 0) return;

	const int trim = 5;

	// Optional: pre-reserve capacity in one go
	const std::size_t max_pairs = n * (n - 1) / 2;
	bambatch_ensure_capacity(R1_batch, max_pairs);
	bambatch_ensure_capacity(R2_batch, max_pairs);

	for (std::size_t i1 = 0; i1 < n; ++i1) {
		for (std::size_t i2 = i1 + 1; i2 < n; ++i2) {

			const std::string& seq1 = frag_seqs[i1];
			const std::string& q1   = frag_quals[i1];

			const std::string& seq2 = frag_seqs[i2];
			const std::string& q2   = frag_quals[i2];

			if (seq2.size() <= static_cast<std::size_t>(trim))
				continue;

			std::string R2_seq  = seq2.substr(trim);
			std::string R2_qual = q2.substr(trim);

			std::string base_id = read.seqHeader + "_" +
								  std::to_string(i1) + "_" +
								  std::to_string(i2);

			std::string R1_id = base_id + " 1";
			std::string R2_id = base_id + " 2";

			bambatch_append_unmapped(
				R1_batch, R1_id, seq1, q1,
				static_cast<uint16_t>(BAM_FPAIRED | BAM_FREAD1)
			);
			bambatch_append_unmapped(
				R2_batch, R2_id, R2_seq, R2_qual,
				static_cast<uint16_t>(BAM_FPAIRED | BAM_FREAD2)
			);
		}
	}
}
