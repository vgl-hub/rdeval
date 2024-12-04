#ifndef READS_H
#define READS_H

#include <iterator>
#include <algorithm>
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
    
    void set(Log* threadLog, uint32_t uId, uint32_t iId, std::string readHeader, std::string* readComment, std::string* read, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, uint32_t readPos, std::string* sequenceQuality, double* avgQuality, std::vector<Tag>* inReadTags = NULL, uint64_t* N = NULL);
    
friend class InReads;

};

class LenVector {

    std::vector<uint8_t> readLens8; // up to 255
    std::vector<uint16_t> readLens16; // up to 65535
    std::vector<uint64_t> readLens64; // max length
    uint64_t curr = 0;
    
public:
    
    void push_back(size_t size) {
        if (size > 255 && size < 65536)
            readLens16.push_back(size);
        else if (size < 256)
            readLens8.push_back(size);
        else
            readLens64.push_back(size);
        
        ++curr;
    }
    
    void insert(const LenVector &vector) {
        readLens8.insert(std::end(readLens8), std::begin(vector.readLens8), std::end(vector.readLens8));
        readLens16.insert(std::end(readLens16), std::begin(vector.readLens16), std::end(vector.readLens16));
        readLens64.insert(std::end(readLens64), std::begin(vector.readLens64), std::end(vector.readLens64));
        curr += vector.readLens8.size() + vector.readLens16.size() + vector.readLens64.size();
    }
    
    uint64_t back() {
        
        uint64_t smallest = 0;
        
        if (readLens8.size())
            smallest = readLens8.back();
        else if (readLens16.size())
            smallest = readLens16.back();
        else if (readLens64.size())
            smallest = readLens64.back();
        
        return smallest;
    }
    
    uint64_t front() {
        
        uint64_t largest = 0;
        
        if (readLens64.size())
            largest = readLens64.front();
        else if (readLens16.size())
            largest = readLens16.back();
        else if (readLens8.size())
            largest = readLens8.back();
        
        return largest;
    }
    
    uint64_t operator[](uint64_t index) {
        
        uint64_t value;
        
        if (index < 0 || index >= (readLens8.size() + readLens16.size() + readLens64.size()))
            throw std::out_of_range("Index out of bounds");
        
        if (index < readLens8.size())
            value = readLens8[index];
        else if (index < (readLens8.size() + readLens16.size()))
            value = readLens16[index-readLens8.size()];
        else
            value = readLens64[index-readLens8.size()-readLens16.size()];
        
        return value;
    }
    
    uint64_t size() {
        return (readLens8.size() + readLens16.size() + readLens64.size());
    }
    
    void sort() {
        std::sort(readLens8.begin(), readLens8.end());
        std::sort(readLens16.begin(), readLens16.end());
        std::sort(readLens64.begin(), readLens64.end());
    }
    
    std::vector<uint64_t> all() { // far from optimal implementation, but needed without iterator
        std::vector<uint64_t> allReadLens;
        allReadLens.reserve(readLens8.size() + readLens16.size() + readLens8.size()); // preallocate memory
        allReadLens.insert(allReadLens.end(), readLens8.begin(), readLens8.end());
        allReadLens.insert(allReadLens.end(), readLens16.begin(), readLens16.end());
        allReadLens.insert(allReadLens.end(), readLens64.begin(), readLens64.end());
        return allReadLens;
    }
    
    std::vector<uint8_t>& getReadLens8() {
        return readLens8;
    }

    std::vector<uint16_t>& getReadLens16() {
        return readLens16;
    }
    
    std::vector<uint64_t>& getReadLens64() {
        return readLens64;
    }
    
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
    
    LenVector readLens;
    
    uint64_t totA=0, totT=0, totC=0, totG=0, totN=0;
    std::vector<std::string> qualities;
    std::vector<double> avgQualities; 
    std::vector<std::vector<int>> qualitiesInts;
    
    OutputStream outputStream;
    uint64_t batchCounter = 1;
    
public:
    
    InReads(UserInputRdeval &userInput, std::string file) : userInput(userInput), outputStream(file) {};
    
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

    void getQualities();

    double getAvgQualities();
    
    void report();

    void printReadLengths();

    void printQualities();

    void printContent();
    
    void evalNstars();
    
    void writeToStream();
    
    void printTableCompressedBinary(std::string outFile);
    
    void readTableCompressedBinary(std::string inFile);
    
};


#endif /* READS_H */
