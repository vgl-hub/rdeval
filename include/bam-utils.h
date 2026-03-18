#ifndef BAMUTILS_H
#define BAMUTILS_H

#include <htslib/sam.h>

#include <string>
#include <vector>
#include <utility>   // std::swap
#include <cstdio>    // fprintf
#include <cstdlib>   // std::abort
#include <cstdint>   // uint16_t, etc

// Create a fresh unmapped BAM record (allocates a new bam1_t).
// NOTE: l_qname passed to bam_set1 must include the trailing '\0' => name.size() + 1.
static inline bam1_t* make_unmapped_bam(
	const std::string& name,
	const std::string& seq,
	const std::string& qual_str,
	uint16_t extra_flags = 0
) {
	bam1_t* b = bam_init1();
	if (!b) {
		std::fprintf(stderr, "Failed to initialize bam1_t\n");
		std::abort();
	}

	const int32_t l_qname = static_cast<int32_t>(name.size() + 1);
	const int32_t l_seq   = static_cast<int32_t>(seq.size());

	const uint16_t flag = static_cast<uint16_t>(BAM_FUNMAP | extra_flags);

	// We pass nullptr for qual here and then fill the numeric quals below.
	// (Some htslib versions allow passing qual_str.c_str() directly; see fill_unmapped_bam)
	if (bam_set1(
			b,
			l_qname,
			name.c_str(),
			flag,
			-1, -1,           // tid, pos
			0,                // mapq
			0,                // n_cigar
			nullptr,          // cigar
			-1, -1,           // mtid, mpos
			0,                // isize
			l_seq,
			seq.c_str(),
			nullptr,          // qual (filled below)
			0                 // l_aux
		) < 0)
	{
		std::fprintf(stderr, "Failed to set BAM data\n");
		std::abort();
	}

	// Fill numeric (0..93) qualities. bam_set1 stores quals as raw Phred, not ASCII.
	uint8_t* qual = bam_get_qual(b);
	if (qual) {
		if (!qual_str.empty()) {
			if (qual_str.size() != seq.size()) {
				std::fprintf(stderr, "Sequence and quality length mismatch in make_unmapped_bam\n");
				std::abort();
			}
			for (size_t i = 0; i < qual_str.size(); ++i) {
				const int q = static_cast<int>(qual_str[i]) - 33;
				qual[i] = static_cast<uint8_t>(q < 0 ? 0 : q);
			}
		} else {
			// If no qualities provided, set to 0 (Phred 0).
			for (int32_t i = 0; i < l_seq; ++i)
				qual[i] = 0;
		}
	}

	return b;
}

// Keep b->data allocated, just mark record empty (pool-friendly).
static inline void bam_recycle_keep_capacity(bam1_t* b) {
	b->l_data = 0;
	b->core = bam1_core_t(); // zero-initialize core
	// Do NOT free b->data; m_data stays as-is.
}

// Fill an existing bam1_t as an unmapped record.
// Uses bam_set1 when available; otherwise falls back to making a temporary record and swapping.
static inline void fill_unmapped_bam(
	bam1_t* b,
	const std::string& name,
	const std::string& seq,
	const std::string& qual_str,
	uint16_t extra_flags = 0
) {
#ifdef bam_set1
	const int32_t l_qname = static_cast<int32_t>(name.size() + 1);
	const int32_t l_seq   = static_cast<int32_t>(seq.size());

	const int ret = bam_set1(
		b,
		l_qname, name.c_str(),
		static_cast<uint16_t>(BAM_FUNMAP | extra_flags),
		-1, -1, 0,              // tid, pos, mapq
		0, nullptr,             // n_cigar, cigar
		-1, -1, 0,              // mtid, mpos, isize
		l_seq,
		seq.c_str(),
		qual_str.empty() ? nullptr : qual_str.c_str(),
		0
	);
	if (ret < 0) {
		std::fprintf(stderr, "bam_set1 failed\n");
		std::abort();
	}

	// If qual_str was nullptr/empty, explicitly set qualities to 0.
	if (qual_str.empty()) {
		uint8_t* q = bam_get_qual(b);
		if (q) {
			for (int32_t i = 0; i < l_seq; ++i)
				q[i] = 0;
		}
	}

#else
	// Fallback: allocate a temporary record using existing helper,
	// then swap contents into b (still avoids repeated malloc/free of b itself).
	bam1_t* tmp = make_unmapped_bam(name, seq, qual_str, extra_flags);

	std::swap(b->core,   tmp->core);
	std::swap(b->l_data, tmp->l_data);
	std::swap(b->m_data, tmp->m_data);
	std::swap(b->data,   tmp->data);

	bam_destroy1(tmp);
#endif
}

// ---- Templated helpers (no dependency on BamBatch typedef) ----
// Requirements on Batch:
//   - batch.reads is std::vector<bam1_t*>
//   - batch.used is an integer count of valid entries

template <typename Batch>
static inline void bambatch_ensure_capacity(Batch& batch, size_t additional) {
	const size_t need = static_cast<size_t>(batch.used) + additional;
	if (need <= batch.reads.size()) return;

	const size_t old = batch.reads.size();
	size_t neu = old ? old : 1024;
	while (neu < need) neu *= 2;

	batch.reads.resize(neu);
	for (size_t i = old; i < neu; ++i)
		batch.reads[i] = bam_init1();
}

template <typename Batch>
static inline void bambatch_append_unmapped(
	Batch& batch,
	const std::string& name,
	const std::string& seq,
	const std::string& qual,
	uint16_t extra_flags
) {
	bambatch_ensure_capacity(batch, 1);

	bam1_t* b = batch.reads[batch.used++];
	bam_recycle_keep_capacity(b);
	fill_unmapped_bam(b, name, seq, qual, extra_flags);
}

#endif /* BAMUTILS_H */
