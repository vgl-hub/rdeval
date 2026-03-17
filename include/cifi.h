#ifndef SCIFI_H
#define SCIFI_H

#include <htslib/sam.h>
#include "bam-utils.h"
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
void chop_read_into_bambatch(const InRead& read, const EnzymeInfo& enz, BamBatch& batch);

void build_pe_bambatches(const InRead& read, const EnzymeInfo& enz, BamBatch& R1_batch, BamBatch& R2_batch);

#endif /* SCIFI_H */
