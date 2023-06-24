#ifndef READS_H
#define READS_H

struct UserInputRdeval : UserInput {

    std::string filter = "none";

};

class InRead : InSegment {

private:
    double avgQuality;
        
public:
    
    void set(Log* threadLog, uint32_t uId, uint32_t iId, std::string readHeader, std::string* readComment, std::string* read, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, uint32_t readPos, std::string* sequenceQuality, double* avgQuality, std::vector<Tag>* inReadTags = NULL, uint64_t* N = NULL);
    
friend class InReads;
    
};

class InReads {
    
    uint32_t batchSize = 10000;
    
    std::vector<Log> logs;
    
    std::shared_ptr<std::istream> stream;
    
    UserInputRdeval userInput;
    std::vector<InRead*> inReads;
    
    //intermediates
    std::string h;
    char* c;
    
    uint32_t seqPos = 0; // to keep track of the original sequence order
    
    std::vector<uint64_t> readNstars    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<uint32_t> readLstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<uint64_t> readLens;
    uint64_t totA =0, totT=0, totC=0, totG=0, totN=0;
    std::vector<std::string> qualities;
    std::vector<double> avgQualities; 
    std::vector<std::vector<int>> qualitiesInts;

    
public:
    
    ~InReads();
    
    void load(UserInputRdeval* userInput);
    
    bool traverseInReads(Sequences* sequence, UserInputRdeval* userInput);
    
    InRead* traverseInRead(Log* threadLog, Sequence* sequence, uint32_t seqPos);
    
    void appendReads(Sequences* readBatch, UserInputRdeval* userInput);
    
    uint64_t getTotReadLen();

    double computeGCcontent();
    
    double computeAvgReadLen();
    
    uint64_t getReadN50();

    uint64_t getSmallestRead();

    uint64_t getLargestRead();

    void getQualities();

    double getAvgQualities();
    
    void report(uint64_t gSize);

    void printReadLengths(char sizeOutType);

    void printQualities(char qualityOut);

    void printContent(char content);
    
    void evalNstars();
    
};


#endif /* READS_H */
