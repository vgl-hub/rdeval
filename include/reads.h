#ifndef READS_H
#define READS_H

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>

#include "len-vector.h"
#include "output.h"

struct UserInputRdeval : UserInput {

    std::string filter = "none";
    
    char sizeOutType = 'u'; //default output from this flag is unsorted sizes
    char qualityOut = 'a'; // average quality per read
    int outSize_flag = 0, quality_flag = 0, content_flag = 0, md5_flag = 0, cmd_flag = 0;
    float ratio = 1.0f;
    int stats_flag = 1; // by default we output the stats
    int32_t randSeed = -1;
};

class InRead : InSegment {

    float avgQuality = 0;
        
public:
    
    void set(Log* threadLog, uint32_t uId, uint32_t iId, std::string readHeader, std::string readComment, std::string read, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, uint32_t readPos, std::string sequenceQuality, float avgQuality, std::vector<Tag>* inReadTags = NULL, uint64_t* N = NULL);
    
friend class InReads;

};

class InReads {
    
    uint32_t batchSize = 1000000; // number of bases processed by a thread
    std::vector<Log> logs;
    std::shared_ptr<std::istream> stream;
    
    UserInputRdeval &userInput;
    std::vector<std::pair<std::vector<InRead*>,uint32_t>> readBatches;
    uint64_t totReads = 0;
    
    //intermediates
    std::string h;
    char* c;
    
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
    bool bam = false;
    
    std::vector<std::pair<std::string,std::string>> md5s;
    
    // dictionaries
    phmap::parallel_flat_hash_set<std::string> includeList;
    phmap::parallel_flat_hash_set<std::string> excludeList;
    
    // filters
    char lSign = '0', qSign = '0', logicalOperator = '0';
    uint64_t l = 0, q = 0;
    
public:
    
    InReads(UserInputRdeval &userInput) : userInput(userInput) {
        
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
    };
    
    void initStream();
    
    void closeStream();
    
    void load();

    void initDictionaries();
    
    void initFilters();
    
    float computeAvgQuality(const std::string &sequenceQuality);
    
    inline bool filterRead(const Sequence &sequence);
    
    inline bool applyFilter(uint64_t size, float avgQuality);
    
    bool traverseInReads(Sequences readBatch);
    
    InRead* traverseInRead(Log* threadLog, Sequence &sequence, uint32_t seqPos);
    
    void appendReads(Sequences readBatch);
    
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
    
    void writeBamHeader();
    
    void closeBam();
};

#endif /* READS_H */
