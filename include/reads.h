#ifndef READS_H
#define READS_H

#include "len-vector.h"
#include "output.h"

struct UserInputRdeval : UserInput {

    std::string filter = "none";
    
    char sizeOutType = 'u'; //default output from this flag is unsorted sizes
    char qualityOut = 'a'; // average quality per read
    char content = 'a'; // default output is to print the normalized ATCGN content for all sequences
    int outSize_flag, quality_flag, content_flag, cmd_flag = 0;

};

class InRead : InSegment {

double avgQuality;
        
public:
    
    void set(Log* threadLog, uint32_t uId, uint32_t iId, std::string readHeader, std::string* readComment, std::string* read, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, uint32_t readPos, std::string* sequenceQuality, double avgQuality, std::vector<Tag>* inReadTags = NULL, uint64_t* N = NULL);
    
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
    
    std::vector<uint64_t> readNstars    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<uint32_t> readLstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    LenVector<float> readLens;
    
    uint64_t totA=0, totT=0, totC=0, totG=0, totN=0;
    std::vector<float> avgQualities; // this is redundant now (already part of the read lens vector), should be removed
    
    OutputStream outputStream;
    bool streamOutput = false;
    uint64_t batchCounter = 1;
    
public:
    
    InReads(UserInputRdeval &userInput, std::string file) : userInput(userInput), outputStream(file) {
        
        const static phmap::flat_hash_map<std::string,int> string_to_case{ // supported read outputs
            {"fasta",1},
            {"fa",1},
            {"fasta.gz",1},
            {"fa.gz",1},
            {"fastq",2},
            {"fq",2},
            {"fastq.gz",2},
            {"fq.gz",2}
        };
        
        if (userInput.outFiles.size()) {
            for (std::string file : userInput.outFiles)
                if (string_to_case.find(getFileExt(file)) != string_to_case.end())
                    streamOutput = true;
        }
    };
    
    void openOutput(std::string file);
    
    void load();
    
    bool traverseInReads(Sequences* sequence);
    
    InRead* traverseInRead(Log* threadLog, Sequence* sequence, uint32_t seqPos);
    
    void appendReads(Sequences* readBatch);
    
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
    
    void printTableCompressedBinary(std::string outFile);
    
    void readTableCompressedBinary(std::string inFile);
    
    void sortReads();
};

#endif /* READS_H */
