#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>

#include <htslib/sam.h>
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

#include "zlib.h"
#include "zstream/zstream_common.hpp"
#include "zstream/ozstream.hpp"
#include "zstream/ozstream_impl.hpp"
#include "output.h"
#include "len-vector.h"

#include "scifi.h"
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
	  free_pool_in(inBuffersN),
	  filled_q_in(inBuffersN),
	  free_pool_out(outBuffersN + outBuffersN * userInput.inputScifi), // two buffers per consumer, double that if scifi (PE output)
	  filled_q_out(outBuffersN + outBuffersN * userInput.inputScifi)
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
		outCount = userInput.scifiCombinations_flag ? numFiles * 2 : numFiles;

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
	
	md5s.reserve(userInput.inFiles.size()); // to avoid invalidating the vector during thread concurrency
	readSummaryBatches.files.resize(userInput.inFiles.size()); // resize to accommodate batches from multiple files
	fileBatches.files.resize(userInput.inFiles.size()); // resize to accommodate batches from multiple files
	
	// Preallocate exactly N buffers
	for (size_t i = 0; i < inBuffersN; ++i) {
		std::unique_ptr<Sequences2> b(new Sequences2);
		free_pool_in.push(std::move(b));
	}
	if (streamOutput) {
		for (size_t i = 0; i < (outBuffersN + outBuffersN * userInput.inputScifi); ++i) {
			auto batch = std::make_unique<BamBatch>();
			batch->reads.reserve(1024);
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
}

void InReads::extractInReads() {
	
	std::string newLine, seqHeader, seqComment, line, bedHeader;
	std::size_t numFiles = userInput.inFiles.size();
	uint32_t batchCounter = 0; // batch number, read position in the batch
	uint64_t processedLength = 0;
	bool sample = userInput.ratio < 1 ? true : false; // read subsampling
	std::unique_ptr<Sequences2> readBatch = free_pool_in.pop();
	
	const static phmap::flat_hash_map<std::string,int> string_to_case{
		{"fasta",1},
		{"fa",1},
		{"fasta.gz",1},
		{"fa.gz",1},
		{"fastq",1},
		{"fq",1},
		{"fastq.gz",1},
		{"fq.gz",1},
		{"bam",2},
		{"cram",2},
		{"rd",3}
	};
    
	while (true) {
		
		uint32_t i = fileCounterStarted.fetch_add(1, std::memory_order_relaxed);
		if (i >= numFiles)
			break;
		
		uint32_t batchN = 0;
        std::string file = userInput.file('r', i);
        std::string ext = getFileExt(file);
        if (ext != "rd") {
            md5s.push_back(std::make_pair(getFileName(file),std::string()));
            threadPool.queueJob([=]{ return computeMd5(file, md5s.back().second); });
        }
        switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
                
            case 1: { // fa*[.gz]
				StreamObj streamObj;
				std::shared_ptr<std::istream> stream = streamObj.openStream(userInput, 'r', i);
				
                if (stream) {
                    
                    switch (stream->peek()) {
                            
                        case '>': {
                            
                            stream->get();
                            
                            while (getline(*stream, newLine)) {
                                
                                size_t spacePos = newLine.find(" ");
                                seqHeader = newLine.substr(0, spacePos);
                                if (spacePos != std::string::npos)
                                    seqComment = newLine.substr(spacePos + 1);
                                else
                                    seqComment.clear();
								std::unique_ptr<std::string> inSequence = std::make_unique<std::string>();
								if (!getline(*stream, *inSequence, '>')) {
									fprintf(stderr, "Record appears truncated (%s). Exiting.\n", seqHeader.c_str());
									exit(EXIT_FAILURE);
								}
								if (sample && newRand() > userInput.ratio)
										continue;
								processedLength += inSequence->size();
								readBatch->add(seqPos++, batchCounter++, std::move(seqHeader), std::move(seqComment), std::move(inSequence));
								
                                if (processedLength > batchSize) {
                                    readBatch->batchN = batchN++;
									readBatch->fileN = i;
                                    lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
									filled_q_in.push(std::move(readBatch));
                                    readBatch = free_pool_in.pop();
                                    processedLength = 0;
									batchCounter = 0;
                                }
                                //lg.verbose("Individual fasta sequence read: " + seqHeader);
                            }
                            break;
                        }
                        case '@': {
                            
                            while (getline(*stream, newLine)) { // file input
                                
                                newLine.erase(0, 1);
                                size_t spacePos = newLine.find(" ");
                                seqHeader = newLine.substr(0, spacePos);
                                if (spacePos != std::string::npos)
                                    seqComment = newLine.substr(spacePos + 1);
                                else
                                    seqComment.clear();
                                
								std::unique_ptr<std::string> inSequence = std::make_unique<std::string>();
								if (!getline(*stream, *inSequence)) {
									fprintf(stderr, "Record appears truncated (%s). Exiting.\n", seqHeader.c_str());
									exit(EXIT_FAILURE);
								}
								if (!getline(*stream, newLine)) {
									fprintf(stderr, "Record appears truncated (%s). Exiting.\n", seqHeader.c_str());
									exit(EXIT_FAILURE);
								}
								std::unique_ptr<std::string> inSequenceQuality = std::make_unique<std::string>();
								if (!getline(*stream, *inSequenceQuality)) {
									fprintf(stderr, "Record appears truncated (%s). Exiting.\n", seqHeader.c_str());
									exit(EXIT_FAILURE);
								}
								if (sample && newRand() > userInput.ratio)
										continue;
								processedLength += inSequence->size();
								readBatch->add(seqPos++, batchCounter++, std::move(seqHeader), std::move(seqComment), std::move(inSequence), std::move(inSequenceQuality));
								
                                if (processedLength > batchSize) {
                                    readBatch->batchN = batchN++;
									readBatch->fileN = i;
                                    lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
									filled_q_in.push(std::move(readBatch));
                                    readBatch = free_pool_in.pop();
                                    processedLength = 0;
									batchCounter = 0;
                                }
                                //lg.verbose("Individual fastq sequence read: " + seqHeader);
                            }
                            break;
                        }
                    }
                    readBatch->batchN = batchN++; // process residual reads
					readBatch->fileN = i;
                    lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
					filled_q_in.push(std::move(readBatch));
					readBatch = free_pool_in.pop();
					processedLength = 0;
					batchCounter = 0;
                }
                break;
            }
            case 2: { // bam, cram
                
                samFile *fp_in = hts_open(file.c_str(),"r"); //open bam file
                bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
                bam1_t *bamdata = bam_init1(); //initialize an alignment
                
				htsThreadPool tpool_read = {NULL, 0};
                tpool_read.pool = hts_tpool_init(userInput.decompression_threads);
                if (tpool_read.pool)
                    hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &tpool_read);
                else
                    lg.verbose("Failed to generate decompression threadpool with " + std::to_string(userInput.decompression_threads) + " threads. Continuing single-threaded");
                
                while(sam_read1(fp_in,bamHdr,bamdata) > 0) {
					
                    uint32_t len = bamdata->core.l_qseq; // length of the read.
                    uint8_t *seq = bam_get_seq(bamdata); // seq string
					std::unique_ptr<std::string> inSequenceQuality;

					std::unique_ptr<std::string> inSequence = std::make_unique<std::string>();
                    inSequence->resize(len);
                    for(uint32_t i=0; i<len; ++i)
                        inSequence->at(i) = seq_nt16_str[bam_seqi(seq,i)]; //gets nucleotide id and converts them into IUPAC id.

					uint8_t* qual = bam_get_qual(bamdata);
					if (qual && (len > 0) && qual[0] != 0xFF) {
						inSequenceQuality = std::make_unique<std::string>(len, '\0');
						for (uint32_t i = 0; i < len; ++i)
							(*inSequenceQuality)[i] = static_cast<char>(qual[i] + 33);
					} else {
						inSequenceQuality = std::make_unique<std::string>(len, '!'); // No per-base qualities; synthesize minimal qualities
					}

					if (sample && newRand() > userInput.ratio)
							continue;
					processedLength += inSequence->size();

					std::string qname(bam_get_qname(bamdata));
					readBatch->add(seqPos++, batchCounter++, std::move(qname), std::string(), std::move(inSequence), std::move(inSequenceQuality));

                    if (processedLength > batchSize) {
                        readBatch->batchN = batchN++;
						readBatch->fileN = i;
                        lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
						filled_q_in.push(std::move(readBatch));
                        readBatch = free_pool_in.pop();
                        processedLength = 0;
						batchCounter = 0;
                    }
                    lg.verbose("Individual fastq sequence read: " + seqHeader);
                }
                readBatch->batchN = batchN++; // process residual reads
				readBatch->fileN = i;
                lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
				filled_q_in.push(std::move(readBatch));
				readBatch = free_pool_in.pop();
				processedLength = 0;
				batchCounter = 0;
                bam_destroy1(bamdata);
                sam_close(fp_in);
                if (tpool_read.pool)
                    hts_tpool_destroy(tpool_read.pool);
                break;
            }
            case 3: { // rd
                    readTableCompressed(userInput.inFiles[i]);
                break;
            }
            default: { // fasta[.gz]
                fprintf(stderr, "cannot recognize input (%s). Must be: fasta, fastq, bam, cram.\n", file.c_str());
                exit(EXIT_FAILURE);
            }
        }
		if (fileCounterCompleted.fetch_add(1, std::memory_order_relaxed) + 1 == numFiles) {
			for (size_t i = 0; i < consumersN; ++i) filled_q_in.push(std::unique_ptr<Sequences2>(nullptr)); // sentinels
		}
    }
}

float InReads::computeAvgQuality(std::string &sequenceQuality) {
    double sumQuality = 0;
    for (char &quality : sequenceQuality)
        sumQuality += pow(10, -(double(quality) - 33)/10);
    return (float) -10 * std::log10(sumQuality/sequenceQuality.size());
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
    
    uint64_t size = sequence->sequence->size();
    float avgQuality = 0;
    if (sequence->sequenceQuality != NULL)
        avgQuality = computeAvgQuality(*sequence->sequenceQuality);
    
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

bool InReads::traverseInReads(Sequences2& readBatchIn)
{
	Log threadLog;
	threadLog.setId(readBatchIn.batchN);
	
	std::unique_ptr<BamBatch> readBatchOut;
	std::unique_ptr<BamBatch> R1_batch;
	std::unique_ptr<BamBatch> R2_batch;

	if (streamOutput) {
		if (userInput.inputScifi && userInput.scifiCombinations_flag) {
			// Two output batches per input batch
			const uint32_t fileIdx = readBatchIn.fileN;
			const uint32_t baseOut = fileIdx * 2;

			R1_batch = free_pool_out.pop();
			R2_batch = free_pool_out.pop();

			R1_batch->reads.clear();
			R2_batch->reads.clear();

			R1_batch->batchN = readBatchIn.batchN;
			R1_batch->fileN  = baseOut;

			R2_batch->batchN = readBatchIn.batchN;
			R2_batch->fileN  = baseOut + 1;
		} else {
			// existing single-output case
			readBatchOut = free_pool_out.pop();
			readBatchOut->reads.clear();
			readBatchOut->fileN  = readBatchIn.fileN;
			readBatchOut->batchN = readBatchIn.batchN;
		}
	}
	EnzymeInfo enz;
	if(userInput.inputScifi)
		enz = get_enzyme(userInput.restrictionEnzyme);
	
	std::vector<InRead> inReadsSummaryBatch;
	uint32_t readN = 0;
	LenVector<float> readLensBatch;

	uint64_t batchA = 0, batchT = 0, batchC = 0, batchG = 0, batchN_ = 0;

	const bool hasInclude = !includeList.empty();
	const bool hasExclude = !excludeList.empty();
	const bool hasFilter  = (userInput.filter != "none");

	for (const auto& sequencePtr : readBatchIn.sequences) {
		Sequence2* sequence = sequencePtr.get();

		if (hasInclude &&
			includeList.find(sequence->header) == includeList.end())
			continue;

		if (hasExclude &&
			excludeList.find(sequence->header) != excludeList.end())
			continue;

		if (hasFilter && filterRead(sequence))
			continue;

		InRead read = traverseInRead(&threadLog,
									 sequence,
									 readBatchIn.batchN + readN++);

		const uint64_t len = read.A + read.C + read.G + read.T + read.N;
		readLensBatch.push_back(std::make_pair(len, read.avgQuality));

		batchA += read.A;
		batchT += read.T;
		batchC += read.C;
		batchG += read.G;
		batchN_ += read.N;

		// STREAMING OUTPUT: convert to BAM record and append to batch
		if (streamOutput && !userInput.inputScifi) {
			// Ensure quality exists
			if (!read.inSequenceQuality) {
				read.inSequenceQuality =
					std::make_unique<std::string>(read.inSequence->size(), '!');
			}

			bam1_t* q = make_unmapped_bam(
				read.seqHeader,              // name
				*read.inSequence,            // seq
				*read.inSequenceQuality,     // qual (ASCII+33)
				0                            // extra_flags (none, just BAM_FUNMAP)
			);

			readBatchOut->reads.push_back(q);
			
		}else if (streamOutput && userInput.inputScifi) {
			
			if (userInput.scifiCombinations_flag) {
				// SCIFI + COMBINATIONS: two *separate* BamBatch outputs

				if (!read.inSequenceQuality) {
					read.inSequenceQuality =
						std::make_unique<std::string>(read.inSequence->size(), '!');
				}
				build_pe_bambatches(read, enz, *R1_batch, *R2_batch);
			}else{
				// SCIFI, NO COMBINATIONS: just chop into fragments in this batch
				if (!read.inSequenceQuality) {
					read.inSequenceQuality =
					std::make_unique<std::string>(read.inSequence->size(), '!');
				}
				
				BamBatch tmp = chop_read_to_bambatch(
													 read,
													 enz,
													 readBatchOut->batchN,
													 readBatchOut->fileN
													 );
				
				// append chopped fragments into our current batch
				readBatchOut->reads.insert(readBatchOut->reads.end(),
										   tmp.reads.begin(), tmp.reads.end());
			}
		}
		if (userInput.content_flag)
			inReadsSummaryBatch.push_back(std::move(read));

		sequence->sequence.reset();
		sequence->sequenceQuality.reset();
	}	
	if (streamOutput) {
		if (userInput.inputScifi && userInput.scifiCombinations_flag) {
				filled_q_out.push(std::move(R1_batch));
				filled_q_out.push(std::move(R2_batch));
			// If they’re empty, they just get destroyed and their BamBatch
			// instances are removed from the pool; the writer will add more
			// via free_pool_out.push after writing other batches.
		} else {
			filled_q_out.push(std::move(readBatchOut));
		}
	}

	{ 	// Update global stats and summaries
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
		totA    += batchA;
		totT    += batchT;
		totC    += batchC;
		totG    += batchG;
		totN    += batchN_;
		totReads += readN;

		threadLog.print();
	}
	writerMutexCondition.notify_one();
	return true;
}

InRead InReads::traverseInRead(Log* threadLog, Sequence2* sequence, uint32_t seqPos) {
	std::vector<std::pair<uint64_t, uint64_t>> bedCoords;
	if (userInput.hc_cutoff != -1) {
		homopolymerCompress(sequence->sequence.get(), bedCoords, userInput.hc_cutoff);
		sequence->sequenceQuality.reset();
	}

	uint64_t A=0, C=0, G=0, T=0, N=0, lowerCount=0;
	for (char base : *sequence->sequence) {
		if (std::islower(static_cast<unsigned char>(base))) ++lowerCount;
		switch (base) {
			case 'A': case 'a': ++A; break;
			case 'C': case 'c': ++C; break;
			case 'G': case 'g': ++G; break;
			case 'T': case 't': ++T; break;
			case 'N': case 'n': case 'X': case 'x': ++N; break;
		}
	}

	float avgQuality = sequence->sequenceQuality
					 ? computeAvgQuality(*sequence->sequenceQuality)
					 : 0.0f;

	InRead inRead(
		threadLog, 0, 0,
		sequence->header, sequence->comment,
		std::move(sequence->sequence),
		std::move(sequence->sequenceQuality),
		seqPos, A, C, G, T, N, lowerCount,
		avgQuality
	);

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
        std::cout<<output("# reads")<<totReads<<std::endl;
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
		outCount = userInput.scifiCombinations_flag ? numFiles * 2 : numFiles;

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
	// Decide output file name
	std::string outName;
	if (splitOutputByFile) { // One/two outputs per input file; typically prefix + per-file suffix
		if (userInput.scifiCombinations_flag) {

			if (outId >= outCount) {
				lg.verbose("openOutputForFile: outId " + std::to_string(outId) +
						   " has no corresponding outFiles entry in scifiCombinations mode");
				std::abort();
			}
			// Determine input file index
			size_t fileIdx = outId / 2;
			if (fileIdx >= numFiles) { // Validate
				lg.verbose("openOutputForFile: computed fileIdx=" + std::to_string(fileIdx) +
						   " invalid (outId=" + std::to_string(outId) +
						   ", total files=" + std::to_string(userInput.inFiles.size()) + ")");
				std::abort();
			}
			std::string baseName = getFileName(userInput.inFiles[fileIdx]);	// Base name of the input file
			std::string ext = getFileExt(userInput.inFiles[fileIdx]);
			std::string baseNameNoExt = stripKnownExt(baseName, ext);
			std::string suffix = ((outId % 2) == 0) ? "_1" : "_2"; // Assign _1 or _2 depending on parity
			outName = userInput.outPrefix + baseNameNoExt + suffix + "." + ext;
			lg.verbose("SciFi combinations mode: outId=" + std::to_string(outId) +
					   " fileIdx=" + std::to_string(fileIdx) +
					   " → " + outName);
		}else{
			if (outId >= userInput.inFiles.size()) {
				lg.verbose("openOutputForFile: outId " + std::to_string(outId) +
						   " has no corresponding outFiles entry");
				std::abort();
			}
			outName = userInput.outPrefix + getFileName(userInput.inFiles[outId]);
		}
	} else { // Single-output mode => everything maps to slot 0
		if (userInput.outFiles.empty()) {
			lg.verbose("openOutputForFile: no output file specified");
			std::abort();
		}
		outName = userInput.outFiles[0];
	}
	// Choose mode based on extension
	std::string ext = getFileExt(outName);

	int fmt_case = string_to_case.count(ext) ? string_to_case.at(ext) : 0;
	htsFile* ofp = nullptr;
	
	switch (fmt_case) {
		case 1:   // fasta[.gz]
		case 2: { // fastq[.gz]
			char mode[4] = "w";
			if (sam_open_mode(mode + 1, outName.c_str(), NULL) < 0) {
				printf("Invalid file name\n");
				exit(EXIT_FAILURE);
			}
			if (!(ofp = sam_open(outName.c_str(), mode))) {
				printf("Could not open %s\n", outName.c_str());
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 3: {  // bam
			ofp = sam_open(outName.c_str(),"wb");
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

	// Duplicate and write BAM/CRAM header
	bam_hdr_t* ohdr = bam_hdr_dup(hdrTemplate);
	if (!ohdr) {
		lg.verbose("Failed to duplicate BAM header");
		std::abort();
	}

	// Attach thread pool if available
	if (tpool_write.pool) {
		hts_set_opt(ofp, HTS_OPT_THREAD_POOL, &tpool_write);
	}

	if (sam_hdr_write(ofp, ohdr) < 0) {
		lg.verbose("Failed to write BAM/CRAM header to " + outName);
		std::abort();
	}

	fps[outId]  = ofp;
	hdrs[outId] = ohdr;

	lg.verbose("Opened output file '" + outName +
			   "' for outId " + std::to_string(outId));
}

void InReads::writeToStream() {
	const size_t numFiles = userInput.inFiles.size();
	size_t outCount = 1;
	if (splitOutputByFile)
		outCount = userInput.scifiCombinations_flag ? numFiles * 2 : numFiles;
	
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
	std::vector<uint64_t> nextBatch(outCount, 0);
	std::vector<std::map<uint64_t, std::unique_ptr<BamBatch>>> pending(outCount);

	for (;;) {
		std::unique_ptr<BamBatch> batch = filled_q_out.pop();

		// nullptr sentinel => no more batches
		if (!batch) break;

		const uint32_t fileId = batch->fileN;
		const uint64_t id     = batch->batchN;

		if (fileId >= outCount) {
			// Defensive: bad file index, recycle (shouldn't happen)
			for (bam1_t* rec : batch->reads)
				bam_destroy1(rec);
			batch->reads.clear();
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

			lg.verbose("Writing " + std::to_string(bptr->reads.size()) +
					   " reads (batch " + std::to_string(bptr->batchN) +
					   ") to outId " + std::to_string(outId) +
					   " (fileId " + std::to_string(fileId) + ")");

			for (bam1_t* rec : bptr->reads) {
				if (sam_write1(ofp, ohdr, rec) < 0) {
					fprintf(stderr, "Error writing BAM record (outId %zu)\n", outId);
					std::abort();
				}
				bam_destroy1(rec);
			}

			bptr->reads.clear();
			free_pool_out.push(std::move(bptr));

			++next;
			it = filePending.find(next);
		}
	}

	// Final flush of any leftover pending batches (should be rare)
	for (size_t fileId = 0; fileId < outCount; ++fileId) {
		auto& filePending = pending[fileId];
		for (auto& kv : filePending) {
			std::unique_ptr<BamBatch>& bptr = kv.second;

			const size_t outId = splitOutputByFile ? fileId : 0;
			openOutputForFile(outId);

			htsFile*   ofp  = fps[outId];
			bam_hdr_t* ohdr = hdrs[outId];

			for (bam1_t* rec : bptr->reads) {
				if (sam_write1(ofp, ohdr, rec) < 0) {
					fprintf(stderr, "Error writing BAM record during final flush (outId %zu)\n", outId);
					std::abort();
				}
				bam_destroy1(rec);
			}

			bptr->reads.clear();
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

void InReads::printTableCompressed(std::string outFile) {
	
    // compute buffer size
    std::vector<std::pair<uint8_t,float>> &readLens8 = readLens.getReadLens8();
    std::vector<std::pair<uint16_t,float>> &readLens16 = readLens.getReadLens16();
    std::vector<std::pair<uint64_t,float>> &readLens64 = readLens.getReadLens64();
    uint64_t len8 = readLens8.size(), len16 = readLens16.size(), len64 = readLens64.size();
    
    uLong sourceLen = sizeof(uint64_t) * (5 + 3) + len8 * sizeof(readLens8[0]) + len16 * sizeof(readLens16[0]) + len64 * sizeof(readLens64[0]); // ACGTN + len8,len16,len64 + vectors
    Bytef *source = new Bytef[sourceLen];
    uLong destLen = compressBound(sourceLen);
    Bytef *dest = new Bytef[destLen];
	std::memset(source, 0, sourceLen);

    unsigned char* ptr = source;

    // write ACGTN counts
    memcpy(ptr, &totA, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &totC, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &totG, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &totT, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &totN, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    
    // write vector lengths
    memcpy(ptr, &len8, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &len16, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &len64, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
 
    // write vectors
    memcpy(ptr, &readLens8[0], len8 * sizeof(readLens8[0]));
    ptr += len8 * sizeof(readLens8[0]);
    memcpy(ptr, &readLens16[0], len16 * sizeof(readLens16[0]));
    ptr += len16 * sizeof(readLens16[0]);
    memcpy(ptr, &readLens64[0], len64 * sizeof(readLens64[0]));
    ptr += len64 * sizeof(readLens64[0]);

	compress2(dest, &destLen, source, sourceLen, Z_BEST_COMPRESSION);
    delete[] source;
    
    std::ofstream ofs(outFile, std::fstream::trunc | std::ios::out | std::ios::binary);
    
    // write md5s
    uint32_t md5sN = md5s.size();
    uint16_t stringSize;
    ofs.write(reinterpret_cast<const char*>(&md5sN), sizeof(uint32_t)); // <-- these could presumably be also gzipped (need to work on the R interface too)
    for (auto md5 : md5s) {
        stringSize = md5.first.size();
        ofs.write(reinterpret_cast<const char*>(&stringSize), sizeof(uint16_t));
        ofs.write(reinterpret_cast<const char*>(md5.first.c_str()), sizeof(char) * stringSize);
        stringSize = md5.second.size();
        ofs.write(reinterpret_cast<const char*>(&stringSize), sizeof(uint16_t));
        ofs.write(reinterpret_cast<const char*>(md5.second.c_str()), sizeof(char) * stringSize);
    }
    ofs.write(reinterpret_cast<const char*>(&sourceLen), sizeof(uint64_t)); // output the decompressed file size
    ofs.write(reinterpret_cast<const char*>(dest), destLen * sizeof(Bytef)); // output compressed data
    ofs.close();
    delete[] dest;
}

void InReads::readTableCompressed(std::string inFile) {
	
	if (userInput.filter != "none" && getFileExt(userInput.file('r', 0)) == "rd")
		filterRecords();
    
    // read
    std::ifstream ifs(inFile, std::ios::binary | std::ios::ate); // compute file size
    std::streamsize fileSize = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    
    uint32_t md5sN;
    uint16_t stringSize;
    ifs.read(reinterpret_cast<char*>(&md5sN), sizeof(uint32_t));
    
    for (uint32_t i = 0; i < md5sN; ++i) {
        std::string filename, md5;
        ifs.read(reinterpret_cast<char*>(&stringSize), sizeof(uint16_t));
        filename.resize(stringSize);
        ifs.read(reinterpret_cast<char*>(&filename[0]), sizeof(char) * stringSize);
        ifs.read(reinterpret_cast<char*>(&stringSize), sizeof(uint16_t));
        md5.resize(stringSize);
        ifs.read(reinterpret_cast<char*>(&md5[0]), sizeof(char) * stringSize);
        md5s.push_back(std::make_pair(filename,md5));
    }
    
    uLongf decompressedSize;
    ifs.read(reinterpret_cast<char*> (&decompressedSize), sizeof(uint64_t)); // read gz-uncompressed size
    Bytef *data = new Bytef[decompressedSize];
    
    uLong compressedSize = fileSize - sizeof(uint64_t); // it seems this should subtract the header too, probably it should be fixed subtracting the md5s
    Bytef *gzData = new Bytef[compressedSize];
    ifs.read(reinterpret_cast<char*> (gzData), compressedSize * sizeof(Bytef)); // read gzipped data
    uncompress(data, &decompressedSize, gzData, compressedSize);
    delete[] gzData;
    
    // read summary statistics
    unsigned char* ptr = data;
    
    uint64_t A, C, G, T, N;
    memcpy(&A, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&C, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&G, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&T, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&N, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    totA += A;
    totC += C;
    totG += G;
    totT += T;
    totN += N;

    uint64_t len8, len16, len64;
    memcpy(&len8, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&len16, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&len64, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);

    // tmp vectors
    LenVector<float> readLensTmp;

    std::vector<std::pair<uint8_t,float>> &readLensTmp8 = readLensTmp.getReadLens8();
    std::vector<std::pair<uint16_t,float>> &readLensTmp16 = readLensTmp.getReadLens16();
    std::vector<std::pair<uint64_t,float>> &readLensTmp64 = readLensTmp.getReadLens64();

    readLensTmp8.resize(len8);
    readLensTmp16.resize(len16);
    readLensTmp64.resize(len64);
    // prefer this to memcpy as std::copy is safer
    std::pair<uint8_t,float>* p8 = reinterpret_cast<std::pair<uint8_t,float>*>(ptr);
    std::copy(p8, p8+len8, readLensTmp8.begin());
    ptr += len8 * sizeof(readLensTmp8[0]);
    std::pair<uint16_t,float>* p16 = reinterpret_cast<std::pair<uint16_t,float>*>(ptr);
    std::copy(p16, p16+len16, readLensTmp16.begin());
    ptr += len16 * sizeof(readLensTmp16[0]);
    std::pair<uint64_t,float>* p64 = reinterpret_cast<std::pair<uint64_t,float>*>(ptr);
    std::copy(p64, p64+len64, readLensTmp64.begin());
    ptr += len64 * sizeof(readLensTmp64[0]);
    
    delete[] data;
    
    // add to vector
    readLens.insert(readLensTmp);

    totReads += len8 + len16 + len64;
}

void InReads::printMd5() {
    for (auto md5 : md5s)
        std::cout<<md5.first<<": "<<md5.second<<std::endl;
}
