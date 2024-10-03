#ifndef READS_H
#define READS_H

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

class InReads {
    
    uint32_t batchSize = 10000;
    
    std::vector<Log> logs;
    
    std::shared_ptr<std::istream> stream;
    
    UserInputRdeval &userInput;
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
    
    InReads(UserInputRdeval &userInput) : userInput(userInput) {};
    
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
    
};


#endif /* READS_H */
