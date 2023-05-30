#ifndef READS_H
#define READS_H

struct UserInputRdeval : UserInput {

    std::string filter = "none";

};

class InRead : InSegment {

private:
    double avgQuality;
        
public:
    
    void set(Log* threadLog, unsigned int uId, unsigned int iId, std::string readHeader, std::string* readComment, std::string* read, unsigned long long int* A, unsigned long long int* C, unsigned long long int* G, unsigned long long int* T, unsigned long long int* lowerCount, unsigned int readPos, std::string* sequenceQuality, double* avgQuality, std::vector<Tag>* inReadTags = NULL, unsigned long long int* N = NULL);
    
friend class InReads;
    
};

class InReads {
    
    unsigned int batchSize = 10000;
    
    std::vector<Log> logs;
    
    std::shared_ptr<std::istream> stream;
    
    UserInputRdeval userInput;
    std::vector<InRead*> inReads;
    
    //intermediates
    std::string h;
    char* c;
    
    unsigned int seqPos = 0; // to keep track of the original sequence order
    
    std::vector<unsigned long long int> readNstars    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> readLstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned long long int> readLens;
    unsigned long long int totA =0, totT=0, totC=0, totG=0, totN=0;
    // std::vector<long double> listA, listC, listT, listG, listN;
    std::vector<std::string> qualities;
    std::vector<double> avgQualities; 
    std::vector<std::vector<int>> qualitiesInts;

    
public:
    
    ~InReads();
    
    void load(UserInputRdeval* userInput);
    
    bool traverseInReads(Sequences* sequence, UserInputRdeval* userInput);
    
    InRead* traverseInRead(Log* threadLog, Sequence* sequence, unsigned int seqPos);
    
    void appendReads(Sequences* readBatch, UserInputRdeval* userInput);
    
    unsigned long long int getTotReadLen();

    int computeGCcontent();
    
    double computeAvgReadLen();
    
    unsigned long long int getReadN50();

    int getSmallestRead();

    int getLargestRead();

    void getQualities();

    double getAvgQualities();
    
    void report(unsigned long long int gSize);

    void printReadLengths(char sizeOutType);

    void printQualities(char qualityOut);

    void printContent(char content);
    
    void evalNstars();

    
};


#endif /* READS_H */
