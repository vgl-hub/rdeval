#ifndef READS_H
#define READS_H

struct UserInputRdeval : UserInput {

    std::string filter;

};

class InRead {

private:
    std::string seqHeader;
    std::string seqComment;
    std::string* inRead = NULL;
    std::string* inSequenceQuality = NULL;
    unsigned long long int A = 0, C = 0, G = 0, T = 0, N = 0, lowerCount = 0;
    double avgQuality;
    unsigned int uId = 0, iId = 0, seqPos = 0; // what are the uId and iId ? 
    std::vector<Tag> tags; // tags refers to ? 

    friend class InReads; // added this like InSegment has to Insequences 

        
public:
    
    ~InRead(); // what is this line doing ? bitwise complement operator? 
    
    void set(Log* threadLog, unsigned int uId, unsigned int iId, std::string readHeader, std::string* readComment, std::string* read, unsigned long long int* A, unsigned long long int* C, unsigned long long int* G, unsigned long long int* T, unsigned long long int* lowerCount, unsigned int readPos, std::string* sequenceQuality, double* avgQuality, std::vector<Tag>* inReadTags = NULL, unsigned long long int* N = NULL);
    
    void setReadHeader(std::string* h);
    
    void setReadComment(std::string c);

    void setInRead(std::string* s);
    
    void setInReadQuality(std::string* q);
    
    void setReadTags(std::vector<Tag>* t);

    void setuId(unsigned int i); // absolute id
    
    void setiId(unsigned int i); // temporary id, internal to scaffold
    
    void setReadPos(unsigned int i); // temporary id, internal to scaffold

    std::string getReadHeader();
    
    std::string getReadComment();
    
    std::vector<Tag> getTags();
    
    std::string getInRead(unsigned int start = 0, unsigned int end = 0);
    
    std::string getInReadQuality(unsigned int start = 0, unsigned int end = 0);

    unsigned int getReadPos();
    
    unsigned long long int getReadLen(unsigned long long int start = 0, unsigned long long int end = 0);
    
    unsigned int getuId(); // absolute id
    
    unsigned int getiId(); // temporary id, internal to scaffold
    
    void setACGT(unsigned long long int* a, unsigned long long int* c, unsigned long long int* g, unsigned long long int* t, unsigned long long int* n = NULL);

    void setAvgQuality(double* avgQuality);

    // void setQualitiesInt(std::vector<int>* readQualitiesInt);
    
    void setLowerCount(unsigned long long int* C);
    
    unsigned long long int getA();
    
    unsigned long long int getC();
    
    unsigned long long int getG();
    
    unsigned long long int getT();

    unsigned long long int getN();
    
    unsigned int getLowerCount(unsigned long long int start = 0, unsigned long long int end = 0);
    
    double computeGCcontent();

    double getAvgQuality(); 
    
    bool trimRead(unsigned int start, unsigned int end);
    
    bool rvcpRead();
    
    bool invertRead();

};

class InReads {
    
    std::vector<Log> logs;
    
    UserInputRdeval userInput;
    std::vector<InRead*> inReads;
    
    //intermediates
    std::string h;
    char* c;
    
    unsigned int seqPos = 0; // to keep track of the original sequence order
    
    std::vector<unsigned long long int> readNstars    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> readLstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned long long int> readLens;
    unsigned long long int totA =0, totT=0, totC=0, totG=0;
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
    
    void evalNstars();

    int meanQual(); 

    int medQual();

    
};


#endif /* READS_H */
