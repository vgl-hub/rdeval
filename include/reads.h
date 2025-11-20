#ifndef READS_H
#define READS_H

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>

#include "len-vector.h"
#include "output.h"
#include "blocking-queue.h"

// -------------------------------------------------------------
// UserInputRdeval
// -------------------------------------------------------------

struct UserInputRdeval : UserInput {
	std::string filter = "none", outPrefix = "";

	char sizeOutType = 'u';
	char qualityOut  = 'a';

	int outSize_flag  = 0;
	int quality_flag  = 0;
	int content_flag  = 0;
	int md5_flag      = 0;
	int cmd_flag      = 0;
	int stats_flag    = 1;
	int cifiCombinations_flag = 0;

	float    ratio       = 1.0f;
	int32_t  randSeed    = -1;
	int      maxThreads  = 8;
	uint16_t parallel_files         = 4;
	uint16_t decompression_threads  = 4;
	uint16_t compression_threads    = 6;
	
	// CiFi
	bool inputCifi = false;
	std::string restrictionEnzyme = "";
};


// -------------------------------------------------------------
// InRead structure
// -------------------------------------------------------------

struct InRead {
	std::string seqHeader;
	std::string seqComment;
	std::unique_ptr<std::string> inSequence;
	std::unique_ptr<std::string> inSequenceQuality;

	uint64_t A=0, C=0, G=0, T=0, N=0, lowerCount=0;
	uint32_t uId=0, iId=0, seqPos=0;

	std::vector<Tag> tags;
	std::vector<std::deque<DBGpath>> variants;
	float avgQuality = 0.0f;

	InRead() = default;

	InRead(const InRead&)            = delete;
	InRead& operator=(const InRead&) = delete;
	InRead(InRead&&) noexcept        = default;
	InRead& operator=(InRead&&) noexcept = default;

	InRead(Log* threadLog,
		   uint32_t uId, uint32_t iId,
		   std::string header, std::string comment,
		   std::unique_ptr<std::string> sequence,
		   std::unique_ptr<std::string> sequenceQuality,
		   uint32_t readPos,
		   uint64_t A, uint64_t C, uint64_t G, uint64_t T,
		   uint64_t N, uint64_t lowerCount,
		   float avgQuality,
		   const std::vector<Tag>* tags = nullptr)
		: seqHeader(std::move(header)),
		  seqComment(std::move(comment)),
		  inSequence(std::move(sequence)),
		  inSequenceQuality(std::move(sequenceQuality)),
		  A(A), C(C), G(G), T(T), N(N), lowerCount(lowerCount),
		  uId(uId), iId(iId), seqPos(readPos), avgQuality(avgQuality)
	{
		(void)threadLog;
		if (tags) this->tags = *tags;
#ifdef DEBUG
		if (threadLog) {
			threadLog->add("Processing read: " + seqHeader +
						   " (uId: " + std::to_string(uId) +
						   ", iId: " + std::to_string(iId) + ")");
		}
#endif
	}

	void dropPayload() {
		inSequence.reset();
		inSequenceQuality.reset();
	}
};


// -------------------------------------------------------------
// ReadBatch and container structures
// -------------------------------------------------------------

template<typename T>
struct ReadBatch {
	std::vector<T> reads;
	uint32_t batchN = 0;    // within-file batch index
	uint32_t fileN  = 0;    // file ID
};

template<typename T>
struct ReadBatches {
	std::vector<std::unique_ptr<ReadBatch<T>>> readBatches;
};

template<typename T>
struct FileBatches {
	std::vector<ReadBatches<T>> files;
};

using BamBatch    = ReadBatch<bam1_t*>;
using InReadBatch = ReadBatch<InRead>;

// -------------------------------------------------------------
// InReads class
// -------------------------------------------------------------

class InReads {

private:
	// Configuration
	uint32_t batchSize = 1'000'000;
	std::atomic<uint32_t> fileCounterStarted{0};
	std::atomic<uint32_t> fileCounterCompleted{0};

	std::vector<Log> logs;

	UserInputRdeval& userInput;

	FileBatches<bam1_t*> fileBatches;
	FileBatches<InRead>  readSummaryBatches;

	uint64_t totReads = 0;
	uint32_t seqPos   = 0;

	std::vector<uint64_t> readNstars{0,0,0,0,0,0,0,0,0,0};
	std::vector<uint32_t> readLstars{0,0,0,0,0,0,0,0,0,0};
	LenVector<float>      readLens;

	uint64_t totA = 0, totT = 0, totC = 0, totG = 0, totN = 0;

	bool streamOutput = false;
	std::thread writer;

	std::mutex              writerMutex;
	std::condition_variable writerMutexCondition;

	// htslib
	bam_hdr_t* hdrTemplate = nullptr;          // template header (e.g. from first input)
	std::vector<htsFile*>   fps;              // output handles: [0] or per-file
	std::vector<bam_hdr_t*> hdrs;             // per-output header (dup of hdrTemplate)

	htsThreadPool tpool_write{nullptr, 0};

	// mode: if true, map input fileN -> output slot fileN
	// if false, everything maps to output slot 0
	bool splitOutputByFile = false;

	// filters and dictionaries
	std::vector<std::pair<std::string,std::string>> md5s;

	phmap::parallel_flat_hash_set<std::string> includeList;
	phmap::parallel_flat_hash_set<std::string> excludeList;

	char lSign='0', qSign='0', logicalOperator='0';
	uint64_t l=0, q=0;

	// MPMC system
	size_t   consumersN = 0;
	uint32_t producersN = 1;

	size_t inBuffersN  = 0;
	size_t outBuffersN = 0;

	// Input queues
	BlockingQueue<std::unique_ptr<Sequences2>> free_pool_in;
	BlockingQueue<std::unique_ptr<Sequences2>> filled_q_in;

	// Output queues (BAM batches)
	BlockingQueue<std::unique_ptr<BamBatch>> free_pool_out;
	BlockingQueue<std::unique_ptr<BamBatch>> filled_q_out;
	
	// CiFi stats
	std::atomic<uint64_t> cifiReadN{0};

public:

	explicit InReads(UserInputRdeval& userInput);

	// Initialization / IO
	void initStream();
	void closeStream();
	void load();

	// Filters / dictionaries
	void initDictionaries();
	void initFilters();
	inline bool filterRead(Sequence2* sequence);
	inline bool applyFilter(uint64_t size, float avgQuality);

	// Processing
	float computeAvgQuality(std::string& sequenceQuality);
	void  filterRecords();
	void  extractInReads();
	bool  traverseInReads(Sequences2& readBatch);
	InRead traverseInRead(Log* threadLog, Sequence2* sequence, uint32_t seqPos);

	// Statistics
	uint64_t getTotReadLen();
	double   computeGCcontent();
	double   computeAvgReadLen();
	uint64_t getReadN50();
	uint64_t getSmallestRead();
	uint64_t getLargestRead();
	double   getAvgQuality();

	void report();
	void printReadLengths();
	void printQualities();
	void printContent();
	void evalNstars();

	// Writing
	void openOutputForFile(size_t outId);
	void writeToStream();

	// Table + md5
	void printTableCompressed(std::string outFile);
	void readTableCompressed(std::string inFile);
	void printMd5();
	void writeHeader();
	void closeBam();
};

#endif /* READS_H */
