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

#include "scifi.h"

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

BamBatch chop_read_to_bambatch(const InRead& read, const EnzymeInfo& enz, uint32_t batchN, uint32_t fileN) {
	if (!read.inSequenceQuality) {
		throw std::runtime_error("InRead::inSequenceQuality is null in chop_read_to_bambatch");
	}

	std::vector<std::string> frag_seqs;
	std::vector<std::string> frag_quals;
	std::vector<int> frag_lengths; // optional, you can drop if not needed

	digest_read(read, enz, frag_seqs, frag_quals, &frag_lengths);

	BamBatch batch;
	batch.batchN = batchN;
	batch.fileN  = fileN;

	batch.reads.reserve(frag_seqs.size());

	for (std::size_t i = 0; i < frag_seqs.size(); ++i) {
		std::string frag_id = read.seqHeader + "_frag" + std::to_string(i);

		bam1_t* b = make_unmapped_bam(
			frag_id,
			frag_seqs[i],
			frag_quals[i],
			0 // no extra flags
		);

		batch.reads.push_back(b);

		// If you want fragment lengths here, record frag_lengths[i] somewhere else.
	}
	return batch;
}

void build_pe_bambatches(const InRead& read, const EnzymeInfo& enz, BamBatch& R1_batch, BamBatch& R2_batch) {
	if (!read.inSequenceQuality) {
		throw std::runtime_error("InRead::inSequenceQuality is null in build_pe_bambatches");
	}

	std::vector<std::string> frag_seqs;
	std::vector<std::string> frag_quals;

	digest_read(read, enz, frag_seqs, frag_quals, nullptr);

	const std::size_t n = frag_seqs.size();
	if (n == 0) return;

	const int trim = 5; // as in original fic2pe
	const std::size_t max_pairs = n * (n - 1) / 2;

	R1_batch.reads.reserve(R1_batch.reads.size() + max_pairs);
	R2_batch.reads.reserve(R2_batch.reads.size() + max_pairs);

	for (std::size_t i1 = 0; i1 < n; ++i1) {
		for (std::size_t i2 = i1 + 1; i2 < n; ++i2) {
			const std::string& seq1 = frag_seqs[i1];
			const std::string& seq2 = frag_seqs[i2];

			const std::string& q1   = frag_quals[i1];
			const std::string& q2   = frag_quals[i2];

			std::string R1_seq  = seq1;
			std::string R1_qual = q1;

			std::string R2_seq;
			std::string R2_qual;

			if (seq2.size() > static_cast<std::size_t>(trim)) {
				R2_seq = seq2.substr(trim);
			} else {
				R2_seq.clear();
			}

			if (q2.size() > static_cast<std::size_t>(trim)) {
				R2_qual = q2.substr(trim);
			} else {
				R2_qual.clear();
			}

			std::string base_id = read.seqHeader + "_" +
								  std::to_string(i1) + "_" +
								  std::to_string(i2);

			std::string R1_id = base_id + " 1";
			std::string R2_id = base_id + " 2";

			bam1_t* b1 = make_unmapped_bam(
				R1_id, R1_seq, R1_qual,
				static_cast<uint16_t>(BAM_FPAIRED | BAM_FREAD1)
			);
			bam1_t* b2 = make_unmapped_bam(
				R2_id, R2_seq, R2_qual,
				static_cast<uint16_t>(BAM_FPAIRED | BAM_FREAD2)
			);

			R1_batch.reads.push_back(b1);
			R2_batch.reads.push_back(b2);
		}
	}
}
