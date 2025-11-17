#ifndef SCIFI_H
#define SCIFI_H

#include <htslib/sam.h>
#include "reads.h"

struct EnzymeInfo {
	std::string name;
	std::string motif;
	int cut_offset;
};

// Get enzyme motif + cut offset
EnzymeInfo get_enzyme(const std::string& name);

// Find restriction enzyme cut positions
std::vector<int> find_cut_positions(const InRead& read, const EnzymeInfo& enz);

void digest_read(const InRead& read, const EnzymeInfo& enz, std::vector<std::string>& frag_seqs, std::vector<std::string>& frag_quals, std::vector<int>* frag_lengths = nullptr);

// Digest sequences and corresponding quality strings
BamBatch chop_read_to_bambatch(const InRead& read, const EnzymeInfo& enz, uint32_t batchN = 0, uint32_t fileN  = 0);

void build_pe_bambatches(const InRead& read, const EnzymeInfo& enz, BamBatch& R1_batch, BamBatch& R2_batch);

static bam1_t* make_unmapped_bam(const std::string& name, const std::string& seq, const std::string& qual_str, uint16_t extra_flags = 0) {
	bam1_t* b = bam_init1();
	if (!b) {
		fprintf(stderr, "Failed to initialize bam1_t\n");
		std::abort();
	}

	const int32_t l_qname = static_cast<int32_t>(name.size());
	const int32_t l_seq   = static_cast<int32_t>(seq.size());

	uint16_t flag = BAM_FUNMAP | extra_flags;

	if (bam_set1(b,
				 l_qname,
				 name.c_str(),
				 flag,
				 -1, -1,
				 0, 0,
				 nullptr,
				 -1, -1,
				 0,
				 l_seq,
				 seq.c_str(),
				 nullptr,
				 0) < 0) {
		fprintf(stderr, "Failed to set BAM data\n");
		std::abort();
	}

	uint8_t* qual = bam_get_qual(b);
	if (qual && !qual_str.empty()) {
		if (qual_str.size() != seq.size()) {
			fprintf(stderr, "Sequence and quality length mismatch in make_unmapped_bam\n");
			std::abort();
		}
		for (size_t i = 0; i < qual_str.size(); ++i) {
			qual[i] = static_cast<uint8_t>(qual_str[i] - 33);
		}
	} else if (qual) {
		for (int32_t i = 0; i < l_seq; ++i)
			qual[i] = 0;
	}

	return b;
}

#endif /* SCIFI_H */
