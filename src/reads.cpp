#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/thread_pool.h>

#include "log.h"
#include "global.h"
#include "bed.h"
#include "struct.h"
#include "gfa-lines.h"
#include "uid-generator.h"
#include "gfa.h"
#include "functions.h" // global functions
#include "stream-obj.h"
#include <iomanip>

#include "output.h"
#include "len-vector.h"

#include "cifi.h"
#include "reads.h"

InReads::InReads(UserInputRdeval& ui)
	: userInput(ui),
	  splitOutputByFile(ui.outPrefix != "" ? true : false),
	  consumersN(static_cast<size_t>(ui.maxThreads)),
	  producersN(
		  ui.outFiles.size() &&
		  getFileExt(ui.outFiles.at(0)) != "rd"
			  ? 1u
			  : std::min<uint32_t>(ui.inFiles.size(), ui.parallel_files)
	  ),
	  inBuffersN((consumersN + producersN) * 2),
	  outBuffersN(consumersN * 4 + 1),
	  free_pool_in(inBuffersN,1),
	  filled_q_in(inBuffersN),
	  free_pool_out(outBuffersN + outBuffersN * userInput.inputCifi), // two buffers per consumer, double that if cifi (PE output)
	  filled_q_out(outBuffersN + outBuffersN * userInput.inputCifi)
{
	static const phmap::flat_hash_map<std::string,int> string_to_case{
		{"fasta",1}, {"fa",1}, {"fasta.gz",1}, {"fa.gz",1},
		{"fastq",2}, {"fq",2}, {"fastq.gz",2},{"fq.gz",2},
		{"bam",3},   {"cram",3}
	};

	if (!userInput.outFiles.empty()) {
		for (const std::string& file : userInput.outFiles) {
			auto ext = getFileExt(file);
			if (string_to_case.find(ext) != string_to_case.end())
				streamOutput = true;
		}
	}
	
	if (splitOutputByFile)
		streamOutput = true;

	const size_t numFiles = userInput.inFiles.size();
	size_t outCount = 1;
	if (splitOutputByFile)
		outCount = userInput.cifiCombinations_flag ? numFiles * 2 : numFiles;

	if (outCount > 0) {
		fps.assign(outCount, nullptr);
		hdrs.assign(outCount, nullptr);
	}

	if (!userInput.inBedInclude.empty() || !userInput.inBedExclude.empty())
		initDictionaries();
	if (userInput.filter != "none")
		initFilters();

	if (userInput.randSeed != -1)
		srandom(userInput.randSeed);
	else
		srandom(time(nullptr));

	if (userInput.maxMem != 0)
		batchSize = userInput.maxMem;
}

float newRand() {
	return (static_cast <float> (random()) / static_cast <float> (RAND_MAX));
}

void InReads::load() {
	uint32_t newMd5s = 0;
	for (std::string file : userInput.inFiles) {
		if (getFileExt(file) != "rd")
			++newMd5s;
	}
	md5s.resize(newMd5s); // to avoid invalidating the vector during thread concurrency
	readSummaryBatches.files.resize(userInput.inFiles.size()); // resize to accommodate batches from multiple files
	fileBatches.files.resize(userInput.inFiles.size()); // resize to accommodate batches from multiple files
	
	if (streamOutput) {
		for (size_t i = 0; i < (outBuffersN + outBuffersN * userInput.inputCifi); ++i) {
			auto batch = std::make_unique<BamBatch>();
			batch->reads.resize(1024);
			for (size_t j = 0; j < batch->reads.size(); ++j)
				batch->reads[j] = bam_init1();
			batch->used = 0;
			free_pool_out.push(std::move(batch));
		}
	}
	lg.verbose("Processing " + std::to_string(userInput.inFiles.size()) + " files");
	lg.verbose("Using " + std::to_string(producersN) + " producers");
	lg.verbose("Using " + std::to_string(consumersN) + " consumers");
	
	std::vector<std::thread> producers;
	for (size_t t = 0; t < producersN; ++t) {
		producers.emplace_back([&]{
				extractInReads();
		});
	}
	
	std::vector<std::thread> consumers;
	for (size_t t = 0; t < consumersN; ++t) {
		consumers.emplace_back([&]{
			for (;;) {
				std::unique_ptr<Sequences2> readBatch = filled_q_in.pop(); // may sleep if empty
				if (!readBatch) { break; }
				traverseInReads(*readBatch);
				readBatch->recycle_keep_capacity(); // keep capacity for reuse
				free_pool_in.push(std::move(readBatch));
			}
		});
	}
	for (auto& th : producers) th.join();
	for (auto& th : consumers) th.join();
	filled_q_out.push(std::unique_ptr<BamBatch>()); // sentinel
	closeStream();
	lg.verbose("free_pool_in created=" + std::to_string(free_pool_in.created()) +
		   " max=" + std::to_string(free_pool_in.max()));
}

void InReads::extractInReads() {

	const std::size_t numFiles = userInput.inFiles.size();
	const bool sample = userInput.ratio < 1;

	// Reusable buffers (per worker thread invocation of extractInReads)
	std::string name, comment, seqBuf, qualBuf, qnameBuf;

	auto split_header = [](const char* s, std::string& header, std::string& comm) {
		const char* sp = std::strchr(s, ' ');
		if (!sp) { header.assign(s); comm.clear(); }
		else { header.assign(s, sp - s); comm.assign(sp + 1); }
	};

	auto push_batch = [&](std::unique_ptr<Sequences2>& readBatch, uint32_t& batchN, uint32_t fileN) {
		readBatch->batchN = batchN++;
		readBatch->fileN  = fileN;
		lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
		filled_q_in.push(std::move(readBatch));
		readBatch = free_pool_in.pop();
	};

	auto enable_threads = [&](htsFile* fp) {
		if (fp && userInput.decompression_threads > 1) {
			// works for BGZF; harmless if unsupported for this handle
			(void)hts_set_threads(fp, userInput.decompression_threads);
		}
	};

	auto read_fa_fq = [&](const std::string& file,
						  std::unique_ptr<Sequences2>& readBatch,
						  uint32_t& batchN,
						  uint32_t fileN) {

		htsFile* fp = hts_open(file.c_str(), "r");
		if (!fp) {
			fprintf(stderr, "Failed to open input (%s)\n", file.c_str());
			exit(EXIT_FAILURE);
		}
		enable_threads(fp);

		kstring_t ks = {0, 0, nullptr};

		// Decide FASTA vs FASTQ from first non-empty line
		int ret = hts_getline(fp, '\n', &ks);
		while (ret >= 0 && ks.l == 0) ret = hts_getline(fp, '\n', &ks);

		if (ret < 0) { // empty file: keep your “residual even if empty” behavior outside
			if (ks.s) free(ks.s);
			hts_close(fp);
			return;
		}

		uint64_t processedLength = 0;

		auto flush_if_needed = [&](size_t add) {
			processedLength += add;
			if (processedLength > batchSize) {
				push_batch(readBatch, batchN, fileN);
				processedLength = 0;
			}
		};

		if (ks.s[0] == '>') {
			// FASTA: header line starts with '>', sequence may be wrapped.
			std::string pendingHeader(ks.s + 1);

			while (true) {
				split_header(pendingHeader.c_str(), name, comment);

				// Decide once per record whether we keep it
				const bool keep = (!sample) || (newRand() <= userInput.ratio);

				seqBuf.clear();
				while (true) {
					ret = hts_getline(fp, '\n', &ks);
					if (ret < 0) { pendingHeader.clear(); break; }          // EOF ends record
					if (ks.l > 0 && ks.s[0] == '>') {                       // next record
						pendingHeader.assign(ks.s + 1);
						break;
					}

					// Only build sequence if we keep this record
					if (keep) seqBuf.append(ks.s, ks.l);
				}

				if (keep) {
					Sequence2& rec = readBatch->next_slot();
					rec.set(std::move(name),
							std::move(comment),
							std::move(seqBuf),
							std::string{}, // no qualities in FASTA
							seqPos++);
					flush_if_needed(rec.sequence.size());
				}

				if (pendingHeader.empty()) break;
			}

		} else if (ks.s[0] == '@') {
			// FASTQ: 4-line records (assumes seq/qual are single-line)
			while (true) {
				if (ks.l == 0 || ks.s[0] != '@') {
					fprintf(stderr, "FASTQ malformed (expected '@') in %s\n", file.c_str());
					exit(EXIT_FAILURE);
				}

				split_header(ks.s + 1, name, comment);

				// Decide once per record whether we keep it
				const bool keep = (!sample) || (newRand() <= userInput.ratio);

				// seq line
				if (hts_getline(fp, '\n', &ks) < 0) {
					fprintf(stderr, "Record appears truncated (%s). Exiting.\n", name.c_str());
					exit(EXIT_FAILURE);
				}
				if (keep) seqBuf.assign(ks.s, ks.l);

				// plus line
				if (hts_getline(fp, '\n', &ks) < 0) {
					fprintf(stderr, "Record appears truncated (%s). Exiting.\n", name.c_str());
					exit(EXIT_FAILURE);
				}

				// qual line
				if (hts_getline(fp, '\n', &ks) < 0) {
					fprintf(stderr, "Record appears truncated (%s). Exiting.\n", name.c_str());
					exit(EXIT_FAILURE);
				}
				if (keep) qualBuf.assign(ks.s, ks.l);

				if (keep) {
					Sequence2& rec = readBatch->next_slot();
					rec.set(std::move(name),
							std::move(comment),
							std::move(seqBuf),
							std::move(qualBuf),
							seqPos++);
					flush_if_needed(rec.sequence.size());
				}

				// next header
				ret = hts_getline(fp, '\n', &ks);
				while (ret >= 0 && ks.l == 0) ret = hts_getline(fp, '\n', &ks);
				if (ret < 0) break;
			}
		} else {
			fprintf(stderr, "Cannot recognize text input (%s). Expected '>' (FASTA) or '@' (FASTQ).\n", file.c_str());
			exit(EXIT_FAILURE);
		}

		if (ks.s) free(ks.s);
		hts_close(fp);
	};

	auto read_bam_cram = [&](const std::string& file,
							 std::unique_ptr<Sequences2>& readBatch,
							 uint32_t& batchN,
							 uint32_t fileN) {

		samFile* fp = hts_open(file.c_str(), "r");
		if (!fp) {
			fprintf(stderr, "Failed to open BAM/CRAM (%s)\n", file.c_str());
			exit(EXIT_FAILURE);
		}
		enable_threads(fp);

		bam_hdr_t* hdr = sam_hdr_read(fp);
		if (!hdr) {
			fprintf(stderr, "Failed to read BAM/CRAM header (%s)\n", file.c_str());
			exit(EXIT_FAILURE);
		}

		bam1_t* b = bam_init1();
		if (!b) {
			fprintf(stderr, "Failed to allocate bam1_t\n");
			exit(EXIT_FAILURE);
		}

		uint64_t processedLength = 0;

		while (sam_read1(fp, hdr, b) > 0) {

			if (sample && newRand() > userInput.ratio)
				continue;

			const uint32_t len = b->core.l_qseq;

			qnameBuf.assign(bam_get_qname(b));

			seqBuf.resize(len);
			char* out = len ? &seqBuf[0] : nullptr;
			uint8_t* bseq = bam_get_seq(b);
			for (uint32_t k = 0; k < len; ++k)
				out[k] = seq_nt16_str[bam_seqi(bseq, k)];

			qualBuf.resize(len);
			char* qout = len ? &qualBuf[0] : nullptr;
			uint8_t* bq = bam_get_qual(b);
			if (bq && len > 0 && bq[0] != 0xFF) {
				for (uint32_t k = 0; k < len; ++k)
					qout[k] = static_cast<char>(bq[k] + 33);
			} else if (len) {
				std::memset(qout, '!', len);
			}

			processedLength += len;

			Sequence2& rec = readBatch->next_slot();
			rec.set(std::move(qnameBuf),
					std::string{},  // comment not stored in BAM path
					std::move(seqBuf),
					std::move(qualBuf),
					seqPos++);

			if (processedLength > batchSize) {
				push_batch(readBatch, batchN, fileN);
				processedLength = 0;
			}
		}

		bam_destroy1(b);
		bam_hdr_destroy(hdr);
		hts_close(fp);
	};

	while (true) {

		uint32_t i = fileCounterStarted.fetch_add(1, std::memory_order_relaxed);
		if (i >= numFiles) break;

		std::unique_ptr<Sequences2> readBatch = free_pool_in.pop();
		uint32_t batchN = 0;

		const std::string file = userInput.file('r', i);
		const std::string ext  = getFileExt(file);

		if (ext != "rd") {
			threadPool.queueJob([this, i, file]() {
				std::string md5;
				computeMd5(file, md5);
				md5s[i].first  = getFileName(file);
				md5s[i].second = md5;
				return true;
			});
		}

		if (ext == "rd") {
			readTableCompressed(userInput.inFiles[i]);
		} else if (ext == "bam" || ext == "cram") {
			read_bam_cram(file, readBatch, batchN, i);
			push_batch(readBatch, batchN, i); // residual (even if empty)
		} else if (ext == "fasta" || ext == "fa" || ext == "fasta.gz" || ext == "fa.gz" ||
				   ext == "fastq" || ext == "fq" || ext == "fastq.gz" || ext == "fq.gz") {
			read_fa_fq(file, readBatch, batchN, i);
			push_batch(readBatch, batchN, i); // residual (even if empty)
		} else {
			fprintf(stderr, "cannot recognize input (%s). Must be: fasta/fastq (optionally .gz), bam, cram, rd.\n",
					file.c_str());
			exit(EXIT_FAILURE);
		}

		if (fileCounterCompleted.fetch_add(1, std::memory_order_relaxed) + 1 == numFiles) {
			for (size_t k = 0; k < consumersN; ++k)
				filled_q_in.push(std::unique_ptr<Sequences2>(nullptr));
		}
	}
}

static inline const double* phredProbLUTd() {
	static double lut[94];
	static bool inited = false;
	if (!inited) {
		for (int q = 0; q < 94; ++q) {
			lut[q] = std::pow(10.0, -double(q) / 10.0);
		}
		inited = true;
	}
	return lut;
}

static inline double avg_q_from_phred33_prob_d(const std::string& qstr) {
	if (qstr.empty()) return 0.0;
	const double* lut = phredProbLUTd();
	double sum = 0.0;
	for (unsigned char c : qstr) {
		int q = int(c) - 33;
		if (q < 0) q = 0;
		if (q > 93) q = 93;
		sum += lut[q];
	}
	const double mean = sum / double(qstr.size());
	return -10.0 * std::log10(mean);
}

float InReads::computeAvgQuality(const std::string& sequenceQuality) {
	return (float)avg_q_from_phred33_prob_d(sequenceQuality);
}

void InReads::initDictionaries() {
	
	std::ifstream includeFile(userInput.inBedInclude);
	std::string key;

	while (includeFile >> key)
		includeList.insert(key);

	includeFile.close();
	
	std::ifstream excludeFile(userInput.inBedExclude);

	while (excludeFile >> key)
		excludeList.insert(key);

	excludeFile.close();
}

void InReads::initFilters() {
	
	bool cannotParse = false;
	
	std::istringstream stream(userInput.filter); // get l,q and logical operator in between
	std::string token1;
	std::string token2;
	std::string token3;
	stream >> token1 >> token2 >> token3;

	if ((userInput.filter.find('l') == std::string::npos && userInput.filter.find('q') == std::string::npos) || // no filter found
		(token1.size() < 3) || // first filter too short
		(token2.size() != 0 && token3.size() == 0) || // second filter missing
		(token2.size() == 0 && token3.size() != 0) || // missing operator
		(token2.size() > 1) // malformatted operator
	)
		cannotParse = true;
	
	if (cannotParse){
		fprintf(stderr, "Could not parse filter. Terminating.\n");
		exit(EXIT_FAILURE);
	}
	
	std::string lFilter, qFilter;
	if (token2.size())
		logicalOperator = token2[0];
	
	if (token1[0] == 'l') {
		lFilter = token1;
		lSign = lFilter[1];
		l = stoi(lFilter.substr(2));
	} else if (token1[0] == 'q') {
		qFilter = token1;
		qSign = qFilter[1];
		q = stoi(qFilter.substr(2));
	}
	
	if (token3.size() != 0) {
		if (token3[0] == 'l') {
			lFilter = token3;
			lSign = lFilter[1];
			l = stoi(lFilter.substr(2));
		}else if (token3[0] == 'q') {
			qFilter = token3;
			qSign = qFilter[1];
			q = stoi(qFilter.substr(2));
		}
	}
}

void InReads::filterRecords() {
	
	uint64_t readCount = readLens.size();
	LenVector<float> readLensFiltered;
	for (uint64_t i = 0; i < readCount; ++i) {
		if (!applyFilter(readLens[i].first, readLens[i].second))
			readLensFiltered.push_back(std::make_pair(readLens[i].first, readLens[i].second));
	}
	if (readLens.size() != readLensFiltered.size()) { // if the read counts after filtering is different these counts are invalidated
		totA = 0;
		totC = 0;
		totG = 0;
		totT = 0;
		totN = 0;
	}
	readLens = readLensFiltered; // overwrite the counts
	totReads = readLens.size();
}

inline bool InReads::filterRead(Sequence2* sequence) {
	
	uint64_t size = sequence->sequence.size();
	float avgQuality = 0.0f;
	if (!sequence->sequenceQuality.empty())
		avgQuality = computeAvgQuality(sequence->sequenceQuality);
	return applyFilter(size, avgQuality);
}

inline bool InReads::applyFilter(uint64_t size, float avgQuality) {
	
	bool lFilter = ((lSign == '0') ||
					((lSign == '>') && (size > l)) ||
					((lSign == '<') && (size < l)) ||
					((lSign == '=') && (size == l))
					);
	
	bool qFilter = ((qSign == '0') ||
					((qSign == '>') && (avgQuality > q)) ||
					((qSign == '<') && (avgQuality < q)) ||
					((qSign == '=') && (avgQuality == q))
					);
	
	if ((logicalOperator == '0' && (lSign != '0' && lFilter)) ||
		(logicalOperator == '0' && (qSign != '0' && qFilter)) ||
		(logicalOperator == '|' && (lFilter || qFilter)) ||
		(logicalOperator == '&' && (lFilter && qFilter))
	   )
		return false;
	return true;
}

bool InReads::traverseInReads(Sequences2& readBatchIn) {
	Log threadLog;
	threadLog.setId(readBatchIn.batchN);

	std::unique_ptr<BamBatch> readBatchOut;
	std::unique_ptr<BamBatch> R1_batch;
	std::unique_ptr<BamBatch> R2_batch;

	// Only need output batches if we're streaming or if CiFi builds internal batches.
	// Keep your original behavior for !streamOutput (alloc once) since you said it was fine.
	if (!streamOutput) {
		if (userInput.inputCifi && userInput.cifiCombinations_flag) {
			R1_batch = std::make_unique<BamBatch>();
			R2_batch = std::make_unique<BamBatch>();
		} else {
			readBatchOut = std::make_unique<BamBatch>();
		}
	}

	// Acquire from pool if streaming
	if (userInput.inputCifi && userInput.cifiCombinations_flag) {
		const uint32_t fileIdx = readBatchIn.fileN;
		const uint32_t baseOut = fileIdx * 2;

		if (streamOutput) {
			R1_batch = free_pool_out.pop();
			R2_batch = free_pool_out.pop();
		}

		R1_batch->used = 0;
		R2_batch->used = 0;

		R1_batch->batchN = readBatchIn.batchN;
		R1_batch->fileN  = baseOut;

		R2_batch->batchN = readBatchIn.batchN;
		R2_batch->fileN  = baseOut + 1;
	} else {
		if (streamOutput)
			readBatchOut = free_pool_out.pop();

		readBatchOut->used  = 0;
		readBatchOut->fileN = readBatchIn.fileN;
		readBatchOut->batchN = readBatchIn.batchN;
	}

	EnzymeInfo enz;
	if (userInput.inputCifi)
		enz = get_enzyme(userInput.restrictionEnzyme);

	std::vector<InRead> inReadsSummaryBatch;
	uint32_t readN = 0;
	LenVector<float> readLensBatch;

	uint64_t batchA = 0, batchT = 0, batchC = 0, batchG = 0, batchN_ = 0;

	const bool hasInclude = !includeList.empty();
	const bool hasExclude = !excludeList.empty();
	const bool hasFilter  = (userInput.filter != "none");
	
	thread_local std::string tmpQual;

	for (uint32_t idx = 0; idx < readBatchIn.used; ++idx) {
		Sequence2& sequence = *readBatchIn.slots[idx];

		if (hasInclude && includeList.find(sequence.header) == includeList.end())
			continue;
		if (hasExclude && excludeList.find(sequence.header) != excludeList.end())
			continue;
		if (hasFilter && filterRead(&sequence))
			continue;

		// Compute stats (and optionally payload) via traverseInRead()
		InRead read = traverseInRead(&threadLog, &sequence, readBatchIn.batchN + readN++);

		const uint64_t len = read.A + read.C + read.G + read.T + read.N;
		readLensBatch.push_back(std::make_pair(len, read.avgQuality));

		batchA  += read.A;
		batchT  += read.T;
		batchC  += read.C;
		batchG  += read.G;
		batchN_ += read.N;

		// ---- STREAMING OUTPUT (non-CiFi): fill BAM directly from Sequence2 ----
		if (streamOutput && !userInput.inputCifi) {
			// Ensure there is a slot in the preallocated bam array
			if (readBatchOut->used >= readBatchOut->reads.size()) {
				// Grow once in a while; still pool-friendly.
				const size_t old = readBatchOut->reads.size();
				const size_t neu = old ? old * 2 : 1024;
				readBatchOut->reads.resize(neu);
				for (size_t j = old; j < neu; ++j)
					readBatchOut->reads[j] = bam_init1();
			}

			bam1_t* q = readBatchOut->reads[readBatchOut->used++];
			bam_recycle_keep_capacity(q);

			if (!sequence.sequenceQuality.empty()) {
				fill_unmapped_bam(q, sequence.header, sequence.sequence, sequence.sequenceQuality, 0);
			} else {
				tmpQual.assign(sequence.sequence.size(), '!');
				fill_unmapped_bam(q, sequence.header, sequence.sequence, tmpQual, 0);
			}
		}
		// ---- CiFi paths (need payload in read) ----
		else if (userInput.inputCifi) {

			// Make sure payload exists for CiFi processing
			if (!read.inSequence) {
				// If traverseInRead was called with needPayload=false by mistake,
				// this makes the behavior safe.
				read.inSequence = std::make_unique<std::string>(sequence.sequence);
			}
			if (!read.inSequenceQuality) {
				if (!sequence.sequenceQuality.empty())
					read.inSequenceQuality = std::make_unique<std::string>(sequence.sequenceQuality);
				else
					read.inSequenceQuality = std::make_unique<std::string>(read.inSequence->size(), '!');
			}

			if (userInput.cifiCombinations_flag) {
				build_pe_bambatches(read, enz, *R1_batch, *R2_batch);
			} else {
				const uint32_t before = readBatchOut->used;
				chop_read_into_bambatch(read, enz, *readBatchOut);
				const uint32_t added  = readBatchOut->used - before;
				cifiReadN.fetch_add(added, std::memory_order_relaxed);
			}
		}

		if (userInput.content_flag)
			inReadsSummaryBatch.push_back(std::move(read));
	}

	// Push output batches
	if (userInput.inputCifi && userInput.cifiCombinations_flag) {
		cifiReadN.fetch_add(R1_batch->used + R2_batch->used);

		if (streamOutput) {
			filled_q_out.push(std::move(R1_batch));
			filled_q_out.push(std::move(R2_batch));
		} else {
			R1_batch->used = 0;
			R2_batch->used = 0;
		}
	} else {
		if (streamOutput)
			filled_q_out.push(std::move(readBatchOut));
		else
			readBatchOut->used = 0;
	}

	{   // Global stats / summaries
		std::unique_lock<std::mutex> lck(writerMutex);

		if (userInput.content_flag && !inReadsSummaryBatch.empty()) {
			if (readBatchIn.fileN >= readSummaryBatches.files.size())
				readSummaryBatches.files.resize(readBatchIn.fileN + 1);

			auto& rbVec = readSummaryBatches.files[readBatchIn.fileN].readBatches;
			auto  summaryBatch = std::make_unique<InReadBatch>();
			summaryBatch->reads  = std::move(inReadsSummaryBatch);
			summaryBatch->batchN = readBatchIn.batchN;
			summaryBatch->fileN  = readBatchIn.fileN;
			rbVec.emplace_back(std::move(summaryBatch));
		}

		readLens.insert(readLensBatch);
		totA     += batchA;
		totT     += batchT;
		totC     += batchC;
		totG     += batchG;
		totN     += batchN_;
		totReads += readN;

		threadLog.print();
	}
	writerMutexCondition.notify_one();
	return true;
}

InRead InReads::traverseInRead(Log* threadLog, Sequence2* sequence, uint32_t seqPos) {
	std::vector<std::pair<uint64_t, uint64_t>> bedCoords;

	if (userInput.hc_cutoff != -1) {
		homopolymerCompress(&sequence->sequence, bedCoords, userInput.hc_cutoff);
		// keep capacity, just clear content
		sequence->sequenceQuality.clear();
	}

	uint64_t A=0, C=0, G=0, T=0, N=0, lowerCount=0;
	for (char base : sequence->sequence) {
		if (std::islower(static_cast<unsigned char>(base))) ++lowerCount;
		switch (base) {
			case 'A': case 'a': ++A; break;
			case 'C': case 'c': ++C; break;
			case 'G': case 'g': ++G; break;
			case 'T': case 't': ++T; break;
			case 'N': case 'n': case 'X': case 'x': ++N; break;
		}
	}

	const float avgQuality =
		(!sequence->sequenceQuality.empty())
			? computeAvgQuality(sequence->sequenceQuality)
			: 0.0f;

	// Decide if this read needs payload in InRead.
	// - Needed for CiFi processing
	// - For content summaries, you currently drop payload when !streamOutput, so payload not needed there.
	// If you later decide content summaries should keep sequence, flip this.
	const bool needPayload = userInput.inputCifi;

	std::unique_ptr<std::string> seqPtr;
	std::unique_ptr<std::string> qualPtr;

	if (needPayload) {
		// yes, this allocates — but ONLY for CiFi (or whatever you later enable)
		seqPtr = std::make_unique<std::string>(sequence->sequence);

		if (!sequence->sequenceQuality.empty())
			qualPtr = std::make_unique<std::string>(sequence->sequenceQuality);
		// else leave nullptr; caller can synthesize if needed
	}

	InRead inRead(
		threadLog, 0, 0,
		sequence->header, sequence->comment,
		std::move(seqPtr),
		std::move(qualPtr),
		seqPos,
		A, C, G, T, N, lowerCount,
		avgQuality
	);

	// Preserve your old behavior: if content summaries + no streaming, drop payload
	if (userInput.content_flag && !streamOutput)
		inRead.dropPayload();

	return inRead;
}

uint64_t InReads::getTotReadLen() {
	
	uint64_t sum = totA + totC + totG + totT + totN;
	if (!sum) {
		uint64_t readCount = readLens.size();
		for (uint64_t i = 0; i < readCount; ++i)
			sum += readLens[i].first;
	}
	return sum;
}

double InReads::computeGCcontent() {

	uint64_t totReadLen = totA + totC + totG + totT;
	double GCcontent = (double) (totG+totC)/totReadLen * 100;
	
	return GCcontent;
}

double InReads::computeAvgReadLen() {
	return (double) getTotReadLen()/totReads;
}

uint64_t InReads::getReadN50() {
	return readNstars[4];
}

void InReads::evalNstars() { // clean up once len-vector iterator is available, still expensive

	uint64_t sum = 0, totLen = getTotReadLen();
	uint8_t N = 1;
	for(unsigned int i = 0; i < readLens.size(); ++i) { // for each length
		sum += readLens[i].first; // increase sum
		while (sum >= ((double) totLen / 10 * N) && N<= 10) { // conditionally add length.at or pos to each N/L* bin
			
			readNstars[N-1] = readLens[i].first;
			readLstars[N-1] = i + 1;
			N += 1;
		}
	}
}

uint64_t InReads::getSmallestRead() {
	return readLens.back();
}

uint64_t InReads::getLargestRead() {
	return readLens.front();
}

double InReads::getAvgQuality(){

	double sumQualities = 0, avgQualitiesSize=readLens.size();

	for (uint64_t i = 0; i < avgQualitiesSize; ++i)
		sumQualities += readLens[i].first * pow(10, -readLens[i].second/10);  // sum the qualities normalized by their read length

	return -10 * std::log10(sumQualities/getTotReadLen());
}

void InReads::report() {

	if (totReads > 0) {
		
		if (userInput.filter != "none" && getFileExt(userInput.file('r', 0)) == "rd")
			filterRecords();
		
		readLens.sort();
		
		std::cout << std::fixed; // disables scientific notation
		std::cout << std::setprecision(2); // 2 decimal points
		
		if (!tabular_flag)
			std::cout<<output("+++Read summary+++")<<"\n";
		std::cout<<output("# reads")<<+totReads<<std::endl;
		std::cout<<output("Total read length")<<getTotReadLen()<<std::endl;
		std::cout<<output("Average read length") << gfa_round(computeAvgReadLen())<<std::endl;
		evalNstars(); // read N* statistics
		std::cout<<output("Read N50")<<getReadN50()<<std::endl;
		std::cout<<output("Smallest read length")<<getSmallestRead()<<std::endl;
		std::cout<<output("Largest read length")<<getLargestRead()<<std::endl;
		std::cout<<output("Coverage")<<gfa_round((double)getTotReadLen()/userInput.gSize)<<std::endl;
		std::cout<<output("GC content %")<<gfa_round(computeGCcontent())<<std::endl;
		std::cout<<output("Base composition (A:C:T:G)")<<totA<<":"<<totC<<":"<<totT<<":"<<totG<<std::endl;
		std::cout<<output("Average per base quality")<<std::abs(getAvgQuality())<<std::endl;
		
		if (userInput.inputCifi) {
			if (!tabular_flag)
				std::cout<<output("+++CiFi summary+++")<<"\n";
			std::cout<<output("# read fragments")<<+cifiReadN<<std::endl;
		}
	}
}

void InReads::printReadLengths() {
	
	std::cout << std::fixed; // disables scientific notation
	std::cout << std::setprecision(2); // 2 decimal points
	
	uint64_t readCount = readLens.size();
	bool noFilter = userInput.filter == "none" ? true : false;
	if(userInput.sizeOutType == 's')
		readLens.sort();
	
	if (userInput.sizeOutType == 'u' || userInput.sizeOutType == 's') {
		
		for (uint64_t i = 0; i < readCount; ++i) {
			if (noFilter || !applyFilter(readLens[i].first, readLens[i].second))
				std::cout << readLens[i].first << "\n";
		}
		
	}else if(userInput.sizeOutType == 'h') {
		
		phmap::parallel_flat_hash_map<uint64_t, uint64_t> hist;
		
		for (uint64_t i = 0; i < readCount; ++i) {
			if (noFilter || !applyFilter(readLens[i].first, readLens[i].second))
				++hist[readLens[i].first];
		}
		std::vector<std::pair<uint64_t, uint64_t>> table(hist.begin(), hist.end()); // converts the hashmap to a table
		std::sort(table.begin(), table.end());
		
		for (auto pair : table)
			std::cout<<pair.first<<"\t"<<pair.second<<"\n";
		
	}else if(userInput.sizeOutType == 'c') {

		phmap::parallel_flat_hash_map<uint64_t, uint64_t> hist;
		
		for (uint64_t i = 0; i < readCount; ++i) {
			if (noFilter || !applyFilter(readLens[i].first, readLens[i].second))
				++hist[readLens[i].first];
		}
		std::vector<std::pair<uint64_t, uint64_t>> table(hist.begin(), hist.end());
		std::sort(table.begin(), table.end());
		uint64_t totReadLen = getTotReadLen(), sum = 0;
		
		for (auto pair : table) {
			std::cout<<+pair.first<<"\t"<<+pair.second<<"\t"<<+pair.first*pair.second<<"\t"<<+(totReadLen - sum)<<"\n";
			sum += pair.first*pair.second;
		}
	}
}

void InReads::printQualities() {
	
	std::cout << std::fixed; // disables scientific notation
	std::cout << std::setprecision(2); // 2 decimal points
	
	uint64_t readCount = readLens.size();
	
	bool noFilter = userInput.filter == "none" ? true : false;
	char qualityOut = userInput.qualityOut;

	if (qualityOut == 'q'){
		for (uint64_t i = 0; i < readCount; ++i) {
			if (noFilter || !applyFilter(readLens[i].first, readLens[i].second))
				std::cout << readLens[i].second << "\n";
		}
	}
	else if (qualityOut == 'a') { // a prints read lengths and qualities
		for (uint64_t i = 0; i < readCount; ++i) {
			if (noFilter || !applyFilter(readLens[i].first, readLens[i].second))
				std::cout << readLens[i].first << "," << readLens[i].second << "\n";
		}
	}
}

void InReads::printContent() {

	// Header
	std::cout << "Header\tComment\tLength\tA\tC\tG\tT\tN\tGC\tAverage Quality\n";

	// Iterate over files
	for (auto& fileRB : readSummaryBatches.files) {
		auto& batches = fileRB.readBatches; // vector<std::unique_ptr<InReadBatch>>

		// Collect raw pointers to sort by batchN
		std::vector<InReadBatch*> sorted;
		sorted.reserve(batches.size());
		for (auto& ptr : batches) {
			if (ptr) sorted.push_back(ptr.get());
		}

		std::sort(sorted.begin(), sorted.end(),
				  [](const InReadBatch* a, const InReadBatch* b) {
					  return a->batchN < b->batchN;
				  });

		// Print reads in sorted batch order
		for (const InReadBatch* batch : sorted) {
			for (const InRead& read : batch->reads) {
				const uint64_t A = read.A;
				const uint64_t C = read.C;
				const uint64_t G = read.G;
				const uint64_t T = read.T;
				const uint64_t N = read.N;
				const uint64_t total = A + C + G + T + N;
				const float gc = total
					? gfa_round(static_cast<float>(G + C) / static_cast<float>(total))
					: 0.0f;

				std::cout << read.seqHeader  << '\t'
						  << read.seqComment << '\t'
						  << total           << '\t'
						  << A << '\t' << C << '\t' << G << '\t' << T << '\t' << N << '\t'
						  << gc              << '\t'
						  << read.avgQuality << '\n';
			}
		}
	}
}

void dump_read(bam1_t* b) {
	printf("->core.tid:(%d)\n", b->core.tid);
	printf("->core.pos:(%ld)\n", static_cast<unsigned long>(b->core.pos));
	printf("->core.bin:(%d)\n", b->core.bin);
	printf("->core.qual:(%d)\n", b->core.qual);
	printf("->core.l_qname:(%d)\n", b->core.l_qname);
	printf("->core.flag:(%d)\n", b->core.flag);
	printf("->core.n_cigar:(%d)\n", b->core.n_cigar);
	printf("->core.l_qseq:(%d)\n", b->core.l_qseq);
	printf("->core.mtid:(%d)\n", b->core.mtid);
	printf("->core.mpos:(%ld)\n", static_cast<unsigned long>(b->core.mpos));
	printf("->core.isize:(%ld)\n", static_cast<unsigned long>(b->core.isize));
	if (b->data) {
		printf("->data:");
		int i;
		for (i = 0; i < b->l_data; ++i) {
			printf("%x ", b->data[i]);
		}
		printf("\n");
	}
	if (b->core.l_qname) {
		printf("qname: %s\n",bam_get_qname(b));
	}
	if (b->core.l_qseq) {
		printf("qseq:");
		int i;
		for (i = 0; i < b->core.l_qseq; ++i) {
			printf("%c", seq_nt16_str[bam_seqi(bam_get_seq(b),i)]);
		}
		printf("\n");
		printf("qual:");
		uint8_t *s = bam_get_qual(b);
		for (i = 0; i < b->core.l_qseq; ++i) {
			printf("%c",s[i] + 33);
		}
		printf("\n");

	}

	if (bam_get_l_aux(b)) {
		uint32_t i = 0;
		uint8_t* aux = bam_get_aux(b);

		while (i < bam_get_l_aux(b)) {
			printf("%.2s:%c:",aux+i,*(aux+i+2));
			i += 2;
			switch (*(aux+i)) {
				case 'Z':
					while (*(aux+1+i) != '\0') { putc(*(aux+1+i), stdout); ++i; }
					break;
			}
			putc('\n',stdout);
			++i;++i;
		}
	}
	printf("\n");
}

void InReads::initStream() {

	if (!streamOutput)
		return;

	writeHeader();  // initializes hdrTemplate only (no file IO)

	// Writer thread will open output files lazily via openOutputForFile()
	writer = std::thread(&InReads::writeToStream, this);
}

void InReads::closeStream() {

	if (writer.joinable())
		writer.join();
}

void InReads::writeHeader() {

	const char init_header[] = "@HD\tVN:1.4\tSO:unknown\n";

	// Destroy any previous template header if present
	if (hdrTemplate) {
		bam_hdr_destroy(hdrTemplate);
		hdrTemplate = nullptr;
	}

	hdrTemplate = bam_hdr_init();
	if (!hdrTemplate) {
		fprintf(stderr, "Failed to allocate BAM header\n");
		std::exit(EXIT_FAILURE);
	}

	hdrTemplate->l_text   = static_cast<int>(std::strlen(init_header));
	hdrTemplate->text     = ::strdup(init_header);  // ownership by hdrTemplate
	hdrTemplate->n_targets = 0;
	hdrTemplate->target_len = nullptr;
	hdrTemplate->target_name = nullptr;

	// The header will be written in openOutputForFile() to each fp[outId]
	// using sam_hdr_write(fps[outId], hdrs[outId]);
}

void InReads::closeBam() {

	if (hdrTemplate) {
		bam_hdr_destroy(hdrTemplate);
		hdrTemplate = nullptr;
	}
}

void InReads::openOutputForFile(size_t outId) {
	const size_t numFiles = userInput.inFiles.size();
	size_t outCount = 1;
	if (splitOutputByFile)
		outCount = userInput.cifiCombinations_flag ? numFiles * 2 : numFiles;

	if (outId >= outCount) {
		lg.verbose("openOutputForFile: invalid outId " + std::to_string(outId));
		return;
	}
	if (fps[outId] != nullptr) // already open
		return;

	const static phmap::flat_hash_map<std::string,int> string_to_case{
		{"fasta",   1},
		{"fa",      1},
		{"fasta.gz",1},
		{"fa.gz",   1},
		{"fastq",   2},
		{"fq",      2},
		{"fastq.gz",2},
		{"fq.gz",   2},
		{"bam",     3},
		{"cram",    4}
	};

	// Decide output file name (UNCHANGED)
	std::string outName;
	if (splitOutputByFile) {
		if (userInput.cifiCombinations_flag) {
			if (outId >= outCount) {
				lg.verbose("openOutputForFile: outId " + std::to_string(outId) +
						   " has no corresponding outFiles entry in cifiCombinations mode");
				std::abort();
			}
			size_t fileIdx = outId / 2;
			if (fileIdx >= numFiles) {
				lg.verbose("openOutputForFile: computed fileIdx=" + std::to_string(fileIdx) +
						   " invalid (outId=" + std::to_string(outId) +
						   ", total files=" + std::to_string(userInput.inFiles.size()) + ")");
				std::abort();
			}
			std::string baseName = getFileName(userInput.inFiles[fileIdx]);
			std::string ext = getFileExt(userInput.inFiles[fileIdx]);
			std::string baseNameNoExt = stripKnownExt(baseName, ext);
			std::string suffix = ((outId % 2) == 0) ? "_1" : "_2";
			outName = userInput.outPrefix + baseNameNoExt + suffix + "." + ext;
			lg.verbose("SciFi combinations mode: outId=" + std::to_string(outId) +
					   " fileIdx=" + std::to_string(fileIdx) +
					   " → " + outName);
		} else {
			if (outId >= userInput.inFiles.size()) {
				lg.verbose("openOutputForFile: outId " + std::to_string(outId) +
						   " has no corresponding outFiles entry");
				std::abort();
			}
			outName = userInput.outPrefix + getFileName(userInput.inFiles[outId]);
		}
	} else {
		if (userInput.outFiles.empty()) {
			lg.verbose("openOutputForFile: no output file specified");
			std::abort();
		}
		outName = userInput.outFiles[0];
	}

	// Choose mode based on extension (UNCHANGED)
	std::string ext = getFileExt(outName);
	int fmt_case = string_to_case.count(ext) ? string_to_case.at(ext) : 0;

	htsFile* ofp = nullptr;
	
	auto strip_bgzf_from_mode = [](char* mode) {
		// Remove 'z' and an optional following digit [0-9]
		// Example: "wfz6" -> "wf", "wFz" -> "wF"
		const size_t n = std::strlen(mode);
		size_t w = 0;
		for (size_t r = 0; r < n; ++r) {
			if (mode[r] == 'z') {
				// skip 'z'
				// if next is a digit, skip it too
				if (r + 1 < n && mode[r + 1] >= '0' && mode[r + 1] <= '9') {
					++r;
				}
				continue;
			}
			mode[w++] = mode[r];
		}
		mode[w] = '\0';
	};

	switch (fmt_case) {

		case 1:   // fasta[.gz]
		case 2: { // fastq[.gz]

			char mode[16] = "w";

			// Let htslib decide format (FASTA vs FASTQ), and it may also add compression flags based on suffix
			// e.g. for ".gz" it often injects 'z' (BGZF) into mode.
			if (sam_open_mode(mode + 1, outName.c_str(), NULL) < 0) {
				printf("Invalid file name\n");
				exit(EXIT_FAILURE);
			}

			if (userInput.bgzip_level >= 0) {
				// BGZF requested: keep/ensure 'z' and optionally set level digit
				int lvl = userInput.bgzip_level;

				// “optional level”: if flag present but no explicit level, you said "use minimum"
				// implement that as level 1 (change to 0 if you truly want no compression)
				if (lvl == 0) lvl = 1;

				if (lvl < 0 || lvl > 9) {
					fprintf(stderr, "Invalid bgzip level %d (must be 0–9)\n", lvl);
					exit(EXIT_FAILURE);
				}

				// Ensure mode has 'z' and a level digit.
				// First strip any existing z/level then append fresh.
				strip_bgzf_from_mode(mode);

				const size_t mlen = std::strlen(mode);
				mode[mlen + 0] = 'z';
				mode[mlen + 1] = char('0' + lvl);
				mode[mlen + 2] = '\0';

			} else {
				// BGZF NOT requested: force plain gzip if ".gz"
				// Remove any bgzf injected by sam_open_mode, then add 'g'
				strip_bgzf_from_mode(mode);

				const size_t mlen = std::strlen(mode);
				mode[mlen + 0] = 'g';
				mode[mlen + 1] = '\0';
			}

			if (!(ofp = sam_open(outName.c_str(), mode))) {
				printf("Could not open %s\n", outName.c_str());
				exit(EXIT_FAILURE);
			}

			break;
		}

		case 3: {  // bam
			ofp = sam_open(outName.c_str(), "wb");
			break;
		}

		case 4: {  // cram
			htsFormat fmt4 = {sequence_data, cram, {3, 1}, gzip, 6, NULL};
			hts_parse_format(&fmt4, "cram,no_ref=1");
			ofp = sam_open_format(outName.c_str(), "wc", &fmt4);
			break;
		}
	}

	if (!ofp) {
		lg.verbose("Failed to open output file: " + outName);
		std::abort();
	}

	if (!hdrTemplate) {
		lg.verbose("openOutputForFile: hdrTemplate is null (no BAM header)");
		std::abort();
	}

	// Duplicate and write header (UNCHANGED)
	bam_hdr_t* ohdr = bam_hdr_dup(hdrTemplate);
	if (!ohdr) {
		lg.verbose("Failed to duplicate BAM header");
		std::abort();
	}

	// Attach thread pool if available (UNCHANGED; works for BGZF too)
	if (tpool_write.pool) {
		hts_set_opt(ofp, HTS_OPT_THREAD_POOL, &tpool_write);
	}

	if (sam_hdr_write(ofp, ohdr) < 0) {
		lg.verbose("Failed to write BAM/CRAM header to " + outName);
		std::abort();
	}

	fps[outId]  = ofp;
	hdrs[outId] = ohdr;

	lg.verbose("Opened output file '" + outName + "' for outId " + std::to_string(outId));
}

void InReads::writeToStream() {
	const size_t numFiles = userInput.inFiles.size();
	
	const size_t numStreams = userInput.cifiCombinations_flag ? numFiles * 2 : numFiles;

	size_t outCount = 1;
	if (splitOutputByFile)
		outCount = numStreams;
	
	if (outCount == 0)
		return;
	// Init threadpool once
	tpool_write = {nullptr, 0};
	tpool_write.pool = hts_tpool_init(userInput.compression_threads);
	if (!tpool_write.pool) {
		lg.verbose("Failed to generate compression threadpool with " +
				   std::to_string(userInput.compression_threads) +
				   " threads. Continuing single-threaded");
	}

	if (fps.size() < outCount)  fps.assign(outCount, nullptr);
	if (hdrs.size() < outCount) hdrs.assign(outCount, nullptr);

	// Per-input-file reorder state
	std::vector<uint64_t> nextBatch(numStreams, 0);
	std::vector<std::map<uint64_t, std::unique_ptr<BamBatch>>> pending(numStreams);

	for (;;) {
		std::unique_ptr<BamBatch> batch = filled_q_out.pop();

		// nullptr sentinel => no more batches
		if (!batch) break;

		const uint32_t fileId = batch->fileN;
		const uint64_t id     = batch->batchN;

		if (fileId >= numStreams) {
			batch->used = 0;
			free_pool_out.push(std::move(batch));
			continue;
		}

		auto& filePending = pending[fileId];
		filePending.emplace(id, std::move(batch));

		auto& next = nextBatch[fileId];
		auto  it   = filePending.find(next);

		while (it != filePending.end()) {
			std::unique_ptr<BamBatch> bptr = std::move(it->second);
			filePending.erase(it);

			// Map input fileId -> output slot
			const size_t outId = splitOutputByFile ? static_cast<size_t>(fileId) : 0;
			openOutputForFile(outId);

			htsFile*   ofp  = fps[outId];
			bam_hdr_t* ohdr = hdrs[outId];

			lg.verbose("Writing " + std::to_string(bptr->used) +
					   " reads (batch " + std::to_string(bptr->batchN) +
					   ") to outId " + std::to_string(outId) +
					   " (fileId " + std::to_string(fileId) + ")");
			
			for (uint32_t k = 0; k < bptr->used; ++k) {
				bam1_t* rec = bptr->reads[k];
				if (sam_write1(ofp, ohdr, rec) < 0) {
					fprintf(stderr, "Error writing BAM record (outId %zu)\n", outId);
					std::abort();
				}
			}
			bptr->used = 0;
			free_pool_out.push(std::move(bptr));

			++next;
			it = filePending.find(next);
		}
	}

	// Final flush of any leftover pending batches (should be rare)
	for (size_t fileId = 0; fileId < numStreams; ++fileId) {
		auto& filePending = pending[fileId];
		for (auto& kv : filePending) {
			std::unique_ptr<BamBatch>& bptr = kv.second;

			const size_t outId = splitOutputByFile ? fileId : 0;
			openOutputForFile(outId);

			htsFile*   ofp  = fps[outId];
			bam_hdr_t* ohdr = hdrs[outId];

			for (uint32_t k = 0; k < bptr->used; ++k) {
				bam1_t* rec = bptr->reads[k];
				if (sam_write1(ofp, ohdr, rec) < 0) {
					fprintf(stderr, "Error writing BAM record during final flush (outId %zu)\n", outId);
					std::abort();
				}
			}
			bptr->used = 0;
			free_pool_out.push(std::move(bptr));
		}
		filePending.clear();
	}

	for (size_t outId = 0; outId < outCount; ++outId) {	// Close all outputs
		if (fps[outId]) {
			sam_close(fps[outId]);
			fps[outId] = nullptr;
		}
		if (hdrs[outId]) {
			bam_hdr_destroy(hdrs[outId]);
			hdrs[outId] = nullptr;
		}
	}

	if (tpool_write.pool) {
		hts_tpool_destroy(tpool_write.pool);
		tpool_write.pool = nullptr;
	}
}

static inline void write_u8(unsigned char*& p, uint8_t v) {
	*p++ = v;
}

static inline void write_u16(unsigned char*& p, uint16_t v) {
	std::memcpy(p, &v, sizeof(v));
	p += sizeof(v);
}

static inline void write_u64(unsigned char*& p, uint64_t v) {
	std::memcpy(p, &v, sizeof(v));
	p += sizeof(v);
}

static inline void write_f32_bits(unsigned char*& p, float v) {
	static_assert(sizeof(float) == 4, "float must be 32-bit");
	uint32_t bits;
	std::memcpy(&bits, &v, sizeof(bits));   // preserve exact bit pattern
	std::memcpy(p, &bits, sizeof(bits));
	p += sizeof(bits);
}

void InReads::printTableCompressed(std::string outFile) {

	auto &readLens8  = readLens.getReadLens8();   // vector<pair<uint8_t,float>>
	auto &readLens16 = readLens.getReadLens16();  // vector<pair<uint16_t,float>>
	auto &readLens64 = readLens.getReadLens64();  // vector<pair<uint64_t,float>>

	const uint64_t len8  = readLens8.size();
	const uint64_t len16 = readLens16.size();
	const uint64_t len64 = readLens64.size();

	// On-disk uncompressed layout:
	// 5x u64 (ACGTN) + 3x u64 (len8,len16,len64)
	// + len8  * (u8  + u32(float bits))
	// + len16 * (u16 + u32(float bits))
	// + len64 * (u64 + u32(float bits))
	const uint64_t sourceLen =
		uint64_t(sizeof(uint64_t) * (5 + 3)) +
		uint64_t(len8)  * (sizeof(uint8_t)  + sizeof(uint32_t)) +
		uint64_t(len16) * (sizeof(uint16_t) + sizeof(uint32_t)) +
		uint64_t(len64) * (sizeof(uint64_t) + sizeof(uint32_t));

	std::vector<Bytef> source(sourceLen);
	unsigned char* ptr = reinterpret_cast<unsigned char*>(source.data());

	// ACGTN counts
	write_u64(ptr, totA);
	write_u64(ptr, totC);
	write_u64(ptr, totG);
	write_u64(ptr, totT);
	write_u64(ptr, totN);

	// vector lengths
	write_u64(ptr, len8);
	write_u64(ptr, len16);
	write_u64(ptr, len64);

	// vectors (padding-free)
	for (const auto& e : readLens8) {
		write_u8(ptr, e.first);
		write_f32_bits(ptr, e.second);
	}
	for (const auto& e : readLens16) {
		write_u16(ptr, e.first);
		write_f32_bits(ptr, e.second);
	}
	for (const auto& e : readLens64) {
		write_u64(ptr, e.first);
		write_f32_bits(ptr, e.second);
	}

	assert(uLong(ptr - reinterpret_cast<unsigned char*>(source.data())) == sourceLen);

	// compress
	uLongf destLen = compressBound(sourceLen);
	std::vector<Bytef> dest(destLen);

	const int zrc = compress2(dest.data(), &destLen, source.data(), sourceLen, Z_BEST_COMPRESSION);
	if (zrc != Z_OK) {
		throw std::runtime_error("compress2() failed with code " + std::to_string(zrc));
	}

	// write output
	std::ofstream ofs(outFile, std::fstream::trunc | std::ios::out | std::ios::binary);
	if (!ofs) throw std::runtime_error("Failed to open output file: " + outFile);

	// md5s (unchanged)
	const uint32_t md5sN = static_cast<uint32_t>(md5s.size());
	ofs.write(reinterpret_cast<const char*>(&md5sN), sizeof(uint32_t));

	uint16_t stringSize;
	for (const auto& md5 : md5s) {
		stringSize = static_cast<uint16_t>(md5.first.size());
		ofs.write(reinterpret_cast<const char*>(&stringSize), sizeof(uint16_t));
		ofs.write(md5.first.data(), stringSize);

		stringSize = static_cast<uint16_t>(md5.second.size());
		ofs.write(reinterpret_cast<const char*>(&stringSize), sizeof(uint16_t));
		ofs.write(md5.second.data(), stringSize);
	}

	// store decompressed size then compressed blob
	ofs.write(reinterpret_cast<const char*>(&sourceLen), sizeof(sourceLen));
	ofs.write(reinterpret_cast<const char*>(dest.data()), destLen);

	if (!ofs) throw std::runtime_error("Write failed (disk full? permission issue?)");
}

static inline uint8_t read_u8(const unsigned char*& p) {
	return *p++;
}

static inline uint16_t read_u16(const unsigned char*& p) {
	uint16_t v;
	std::memcpy(&v, p, sizeof(v));
	p += sizeof(v);
	return v;
}

static inline uint64_t read_u64(const unsigned char*& p) {
	uint64_t v;
	std::memcpy(&v, p, sizeof(v));
	p += sizeof(v);
	return v;
}

static inline float read_f32_bits(const unsigned char*& p) {
	static_assert(sizeof(float) == 4, "float must be 32-bit");
	uint32_t bits;
	std::memcpy(&bits, p, sizeof(bits));
	p += sizeof(bits);
	float v;
	std::memcpy(&v, &bits, sizeof(v));
	return v;
}

void InReads::readTableCompressed(std::string inFile) {

	if (userInput.filter != "none" && getFileExt(userInput.file('r', 0)) == "rd")
		filterRecords();

	// open and get file size
	std::ifstream ifs(inFile, std::ios::binary | std::ios::ate);
	if (!ifs) {
		throw std::runtime_error("Failed to open input file: " + inFile);
	}

	const std::streamsize fileSize = ifs.tellg();
	if (fileSize < 0) {
		throw std::runtime_error("tellg() failed for: " + inFile);
	}
	ifs.seekg(0, std::ios::beg);

	// ---- read md5 header and track consumed bytes ----
	std::streamsize headerBytes = 0;

	uint32_t md5sN = 0;
	ifs.read(reinterpret_cast<char*>(&md5sN), sizeof(md5sN));
	if (!ifs) {
		throw std::runtime_error("Failed to read md5sN from: " + inFile);
	}
	headerBytes += static_cast<std::streamsize>(sizeof(md5sN));

	for (uint32_t i = 0; i < md5sN; ++i) {
		uint16_t stringSize = 0;
		std::string filename, md5;

		ifs.read(reinterpret_cast<char*>(&stringSize), sizeof(stringSize));
		if (!ifs) {
			throw std::runtime_error("Failed to read filename length from: " + inFile);
		}
		headerBytes += static_cast<std::streamsize>(sizeof(stringSize));

		filename.resize(static_cast<size_t>(stringSize));
		if (stringSize > 0) {
			ifs.read(&filename[0], static_cast<std::streamsize>(stringSize));
			if (!ifs) {
				throw std::runtime_error("Failed to read filename from: " + inFile);
			}
		}
		headerBytes += static_cast<std::streamsize>(stringSize);

		ifs.read(reinterpret_cast<char*>(&stringSize), sizeof(stringSize));
		if (!ifs) {
			throw std::runtime_error("Failed to read md5 length from: " + inFile);
		}
		headerBytes += static_cast<std::streamsize>(sizeof(stringSize));

		md5.resize(static_cast<size_t>(stringSize));
		if (stringSize > 0) {
			ifs.read(&md5[0], static_cast<std::streamsize>(stringSize));
			if (!ifs) {
				throw std::runtime_error("Failed to read md5 from: " + inFile);
			}
		}
		headerBytes += static_cast<std::streamsize>(stringSize);

		md5s.emplace_back(std::move(filename), std::move(md5));
	}

	// ---- read decompressed size from file format ----
	uint64_t decompressedSize = 0;
	ifs.read(reinterpret_cast<char*>(&decompressedSize), sizeof(decompressedSize));
	if (!ifs) {
		throw std::runtime_error("Failed to read decompressed size from: " + inFile);
	}
	headerBytes += static_cast<std::streamsize>(sizeof(decompressedSize));

	if (decompressedSize == 0) {
		throw std::runtime_error("Invalid decompressed size (0) in: " + inFile);
	}
	if (headerBytes > fileSize) {
		throw std::runtime_error("Corrupt file (header larger than file) in: " + inFile);
	}

	// ---- compressed payload size from file ----
	const std::streamsize compressedSize_ss = fileSize - headerBytes;
	if (compressedSize_ss <= 0) {
		throw std::runtime_error("Invalid compressed payload size in: " + inFile);
	}
	const uint64_t compressedSize64 = static_cast<uint64_t>(compressedSize_ss);

	// ---- validate sizes for allocation ----
	if (decompressedSize > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
		throw std::runtime_error("Decompressed payload too large to allocate in: " + inFile);
	}
	if (compressedSize64 > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
		throw std::runtime_error("Compressed payload too large to allocate in: " + inFile);
	}

	const size_t decompressedSize_sz = static_cast<size_t>(decompressedSize);
	const size_t compressedSize_sz = static_cast<size_t>(compressedSize64);

	// ---- validate sizes for zlib API ----
	if (decompressedSize > static_cast<uint64_t>(std::numeric_limits<uLongf>::max())) {
		throw std::runtime_error("Decompressed payload too large for zlib API in: " + inFile);
	}
	if (compressedSize64 > static_cast<uint64_t>(std::numeric_limits<uLong>::max())) {
		throw std::runtime_error("Compressed payload too large for zlib API in: " + inFile);
	}

	const uLongf decompressedSize_z = static_cast<uLongf>(decompressedSize);
	const uLong compressedSize_z = static_cast<uLong>(compressedSize64);

	// ---- read compressed payload ----
	std::vector<Bytef> compressed(compressedSize_sz);
	ifs.read(reinterpret_cast<char*>(compressed.data()), compressedSize_ss);
	if (!ifs) {
		throw std::runtime_error("Failed to read compressed payload from: " + inFile);
	}

	// ---- decompress ----
	std::vector<Bytef> data(decompressedSize_sz);
	uLongf outSize = decompressedSize_z;

	const int zrc = uncompress(
		data.data(), &outSize,
		compressed.data(), compressedSize_z
	);

	if (zrc != Z_OK) {
		throw std::runtime_error("uncompress() failed with code " + std::to_string(zrc) + " for: " + inFile);
	}
	if (outSize != decompressedSize_z) {
		throw std::runtime_error("uncompress() size mismatch for: " + inFile);
	}

	// ---- parse decompressed buffer ----
	const unsigned char* ptr = reinterpret_cast<const unsigned char*>(data.data());
	const unsigned char* end = ptr + data.size();

	auto need = [&](size_t n) {
		if (static_cast<size_t>(end - ptr) < n) {
			throw std::runtime_error("Corrupt decompressed payload (truncated) in: " + inFile);
		}
	};

	// ACGTN
	need(sizeof(uint64_t) * 5);
	const uint64_t A = read_u64(ptr);
	const uint64_t C = read_u64(ptr);
	const uint64_t G = read_u64(ptr);
	const uint64_t T = read_u64(ptr);
	const uint64_t N = read_u64(ptr);

	totA += A;
	totC += C;
	totG += G;
	totT += T;
	totN += N;

	// counts
	need(sizeof(uint64_t) * 3);
	const uint64_t len8  = read_u64(ptr);
	const uint64_t len16 = read_u64(ptr);
	const uint64_t len64 = read_u64(ptr);

	// validate expected payload size exactly
	const uint64_t expectedPayloadSize =
		uint64_t(sizeof(uint64_t) * (5 + 3)) +
		len8  * uint64_t(sizeof(uint8_t)  + sizeof(uint32_t)) +
		len16 * uint64_t(sizeof(uint16_t) + sizeof(uint32_t)) +
		len64 * uint64_t(sizeof(uint64_t) + sizeof(uint32_t));

	if (expectedPayloadSize != decompressedSize) {
		throw std::runtime_error(
			"Decompressed payload size does not match embedded record counts in: " + inFile
		);
	}

	// tmp vectors
	LenVector<float> readLensTmp;
	auto& readLensTmp8  = readLensTmp.getReadLens8();
	auto& readLensTmp16 = readLensTmp.getReadLens16();
	auto& readLensTmp64 = readLensTmp.getReadLens64();

	if (len8 > static_cast<uint64_t>(std::numeric_limits<size_t>::max()) ||
		len16 > static_cast<uint64_t>(std::numeric_limits<size_t>::max()) ||
		len64 > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
		throw std::runtime_error("Record count too large to allocate in: " + inFile);
	}

	readLensTmp8.reserve(static_cast<size_t>(len8));
	readLensTmp16.reserve(static_cast<size_t>(len16));
	readLensTmp64.reserve(static_cast<size_t>(len64));

	// entries: (u8 + f32 bits)
	for (uint64_t i = 0; i < len8; ++i) {
		need(sizeof(uint8_t) + sizeof(uint32_t));
		const uint8_t l = read_u8(ptr);
		const float   w = read_f32_bits(ptr);
		readLensTmp8.emplace_back(l, w);
	}

	// entries: (u16 + f32 bits)
	for (uint64_t i = 0; i < len16; ++i) {
		need(sizeof(uint16_t) + sizeof(uint32_t));
		const uint16_t l = read_u16(ptr);
		const float    w = read_f32_bits(ptr);
		readLensTmp16.emplace_back(l, w);
	}

	// entries: (u64 + f32 bits)
	for (uint64_t i = 0; i < len64; ++i) {
		need(sizeof(uint64_t) + sizeof(uint32_t));
		const uint64_t l = read_u64(ptr);
		const float    w = read_f32_bits(ptr);
		readLensTmp64.emplace_back(l, w);
	}

	// enforce exact consumption
	if (ptr != end) {
		throw std::runtime_error("Extra bytes in decompressed payload in: " + inFile);
	}

	// merge
	readLens.insert(readLensTmp);
	totReads += len8 + len16 + len64;
}


void InReads::printMd5() {
	for (auto md5 : md5s)
		std::cout<<md5.first<<": "<<md5.second<<std::endl;
}
