#ifndef READS_H
#define READS_H

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>

#include "len-vector.h"
#include "output.h"
#include "blocking-queue.h"

struct UserInputRdeval : UserInput {

    std::string filter = "none";
    
    char sizeOutType = 'u'; //default output from this flag is unsorted sizes
    char qualityOut = 'a'; // average quality per read
    int outSize_flag = 0, quality_flag = 0, content_flag = 0, md5_flag = 0, cmd_flag = 0;
    float ratio = 1.0f;
    int stats_flag = 1; // by default we output the stats
    int32_t randSeed = -1;
	int maxThreads = 8;
    uint8_t decompression_threads = 4, compression_threads = 6;
};

struct InRead {
	std::string seqHeader;
	std::string seqComment;
	std::unique_ptr<std::string> inSequence;
	std::unique_ptr<std::string> inSequenceQuality;
	uint64_t A = 0, C = 0, G = 0, T = 0, N = 0, lowerCount = 0;
	unsigned int uId = 0, iId = 0, seqPos = 0;
	std::vector<Tag> tags;
	std::vector<std::deque<DBGpath>> variants;
	float avgQuality = 0.0f;
	
	InRead() = default;

	// COPY: disable
	InRead(const InRead&) = delete;
	InRead& operator=(const InRead&) = delete;

	// MOVE: enable
	InRead(InRead&&) noexcept = default;
	InRead& operator=(InRead&&) noexcept = default;

	// main constructor
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
		if (tags) this->tags = *tags;
#ifdef DEBUG
		if (threadLog)
			threadLog->add("Processing read: " + seqHeader +
						   " (uId: " + std::to_string(uId) +
						   ", iId: " + std::to_string(iId) + ")");
#endif
	}

	void dropPayload() {
		inSequence.reset();
		inSequenceQuality.reset();
	}
};

class InReads {
	
	uint32_t batchSize = 1000000; // number of bases processed by a thread
	std::atomic<uint32_t> fileCounterStarted{0}, fileCounterCompleted{0}; // next file to read
    std::vector<Log> logs;
    
    UserInputRdeval &userInput;
    std::vector<std::pair<std::vector<bam1_t*>,uint32_t>> readBatches;
	std::vector<std::vector<std::pair<std::vector<InRead>,uint32_t>>> readSummaryBatches; // could be avoided in the future
    uint64_t totReads = 0;
    
    uint32_t seqPos = 0; // to keep track of the original sequence order
    
    std::vector<uint64_t> readNstars{0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<uint32_t> readLstars{0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    LenVector<float> readLens;
    uint64_t totA=0, totT=0, totC=0, totG=0, totN=0;
    
    bool streamOutput = false;
    std::thread writer;
    std::condition_variable writerMutexCondition;
    
    htsFile *fp; // htslib file pointer
    bam_hdr_t *hdr; // htslib sam header pointer
    htsThreadPool tpool_write; // htslib threadpool pointer
    bool bam = false;
    
    std::vector<std::pair<std::string,std::string>> md5s;
    
    // dictionaries
    phmap::parallel_flat_hash_set<std::string> includeList;
    phmap::parallel_flat_hash_set<std::string> excludeList;
    
    // filters
    char lSign = '0', qSign = '0', logicalOperator = '0';
    uint64_t l = 0, q = 0;
	
	// MPMC
	size_t NUM_BUFFERS, QCAP;
	BlockingQueue<std::unique_ptr<Sequences2>> free_pool, filled_q;
	size_t N_CONS = userInput.maxThreads;
    
public:
    
    InReads(UserInputRdeval &userInput) : userInput(userInput), NUM_BUFFERS(userInput.maxThreads*2), QCAP(NUM_BUFFERS), free_pool(QCAP), filled_q(QCAP) {
        
        const static phmap::flat_hash_map<std::string,int> string_to_case{ // supported read outputs
            {"fasta",1},
            {"fa",1},
            {"fasta.gz",1},
            {"fa.gz",1},
            {"fastq",2},
            {"fq",2},
            {"fastq.gz",2},
            {"fq.gz",2},
            {"bam",3},
            {"cram",3}
        };
        
        if (userInput.outFiles.size()) {
            for (std::string file : userInput.outFiles) {
                if (string_to_case.find(getFileExt(file)) != string_to_case.end())
                    streamOutput = true;
                
                if (getFileExt(file) == "bam" || getFileExt(file) == "cram")
                    bam = true;
            }
        }
        if(userInput.inBedInclude != "" || userInput.inBedExclude != "")
            initDictionaries();
        if(userInput.filter != "none")
            initFilters();
        
        
        if (userInput.randSeed != -1)
            srandom(userInput.randSeed);
        else
            srandom(time(nullptr));
		
		if (userInput.maxMem != 0)
			batchSize = userInput.maxMem;
    };
    
    void initStream();
    
    void closeStream();
    
    void load();

    void initDictionaries();
    
    void initFilters();
    
    float computeAvgQuality(std::string &sequenceQuality);
	
	void filterRecords();
    
    inline bool filterRead(Sequence2* sequence);
    
    inline bool applyFilter(uint64_t size, float avgQuality);
	
	void extractInReads();
    
    bool traverseInReads(Sequences2 &readBatch);
    
    InRead traverseInRead(Log* threadLog, Sequence2* sequence, uint32_t seqPos);
    
    uint64_t getTotReadLen();

    double computeGCcontent();
    
    double computeAvgReadLen();
    
    uint64_t getReadN50();

    uint64_t getSmallestRead();

    uint64_t getLargestRead();

    double getAvgQuality();
    
    void report();

    void printReadLengths();

    void printQualities();

    void printContent();
    
    void evalNstars();
    
    void writeToStream();
    
    void printTableCompressed(std::string outFile);
    
    void readTableCompressed(std::string inFile);
    
    void printMd5();
    
    void writeHeader();
    
    void closeBam();
};

#endif /* READS_H */
