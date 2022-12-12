#include <stdlib.h>
#include <string.h>

#include <istream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>

#include "log.h"
#include "global.h"

#include "bed.h"
#include "struct.h"
#include "gfa-lines.h"
#include "functions.h" // global functions

#include "stream-obj.h"

#include "reads.h"

InReads::~InReads()
{

    for (InRead* p : inReads)
        delete p;

}

void InReads::load(UserInputRdeval* userInput) { 

    unsigned int batchSize = 100;
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    std::shared_ptr<std::istream> stream;
    
    unsigned int numFiles = userInput->iReadFileArg.size();
    
    lg.verbose("Processing " + std::to_string(numFiles) + " files");
    
    for (unsigned int i = 0; i < numFiles; i++) {
        
        StreamObj streamObj;

        stream = streamObj.openStream(*userInput, 'r', &i);

        Sequences* readBatch = new Sequences;

        if (stream) {

            switch (stream->peek()) {

                case '>': {

                    stream->get();

                    while (getline(*stream, newLine)) {
                        
                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }

                        std::string* inSequence = new std::string;

                        getline(*stream, *inSequence, '>');

                        readBatch->sequences.push_back(new Sequence {seqHeader, seqComment, inSequence});
                        seqPos++;

                        if (seqPos % batchSize == 0) {

                            readBatch->batchN = seqPos/batchSize;
                            
                            lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));

                            appendReads(readBatch, userInput);

                            readBatch = new Sequences;

                        }

                        lg.verbose("Individual fasta sequence read: " + seqHeader);

                    }

                    break;
                }
                case '@': {

                    while (getline(*stream, newLine)) { // file input

                        newLine.erase(0, 1);

                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }

                        std::string* inSequence = new std::string;
                        getline(*stream, *inSequence);

                        getline(*stream, newLine);

                        std::string* inSequenceQuality = new std::string;
                        getline(*stream, *inSequenceQuality);

                        readBatch->sequences.push_back(new Sequence {seqHeader, seqComment, inSequence, inSequenceQuality});
                        seqPos++;

                        if (seqPos % batchSize == 0) {

                            readBatch->batchN = seqPos/batchSize;
                            
                            lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));

                            appendReads(readBatch, userInput);

                            readBatch = new Sequences;

                        }

                        lg.verbose("Individual fastq sequence read: " + seqHeader);

                    }

                    break;

                }

            }
            
            readBatch->batchN = seqPos/batchSize + 1;
                
            lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));

            appendReads(readBatch, userInput);

        }
        
    }

}

void InReads::appendReads(Sequences* readBatch, UserInputRdeval* userInput) { // read a collection of reads
    
    threadPool.queueJob([=]{ return traverseInReads(readBatch, userInput); });
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
    lck.unlock();
    
}

bool InReads::traverseInReads(Sequences* readBatch, UserInputRdeval* userInput) { // traverse the read

    Log threadLog;
    
    threadLog.setId(readBatch->batchN);
    
    std::vector<InRead*> inReadsBatch;
    
    unsigned int readN = 0; 
    std::vector<unsigned long long int> readLensBatch;
    InRead* read;
    unsigned long long int batchA = 0, batchT=0, batchC=0, batchG=0, batchN =0;
    std::vector<double> batchAvgQualities;
    unsigned int filterInt;
    if (!(userInput->filter == "none")) {
        filterInt = stoi(userInput->filter.substr(1));
    }


    // perhaps I need to hard code a small function because the filter is currently in the for-loop and basically requires re-parsing the filter operand each time

    for (Sequence* sequence : readBatch->sequences) {

        if ((userInput->filter[0] == '>') && (sequence->sequence->size()< filterInt)){ 
            // std::cout << "Sequence shorter than filter length." << "\n"; // less compute time to assign to a variable or call like this? // would like to output sequence name here 
            continue;
        }
        else if ((userInput->filter[0] == '<') && (sequence->sequence->size() > filterInt)) {
            // std::cout << "Sequence longer than filter length." << "\n"; 
            continue;
        }
        else if ((userInput->filter[0] == '=') && (sequence->sequence->size() != filterInt)) {
            // std::cout << "Sequence does not equal filter length." << "\n";
            continue;
        }
        
        read = traverseInRead(&threadLog, sequence, readBatch->batchN+readN++);

        readLensBatch.push_back(read->getReadLen());
        batchA += read->getA();
        batchT += read->getT();
        batchC += read->getC();
        batchG += read->getG();
        batchN += read->getN();

        if (read->inSequenceQuality != NULL) 
            batchAvgQualities.push_back(read->getAvgQuality());
        
        
        inReadsBatch.push_back(read);
        
    }
    
    delete readBatch;
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    inReads.insert(std::end(inReads), std::begin(inReadsBatch), std::end(inReadsBatch));
    readLens.insert(std::end(readLens), std::begin(readLensBatch), std::end(readLensBatch));
    avgQualities.insert(std::end(avgQualities), std::begin(batchAvgQualities), std::end(batchAvgQualities));


    totA+=batchA;
    totT+=batchT;
    totC+=batchC;
    totG+=batchG;
    totN+=batchN;

    logs.push_back(threadLog);
    
    lck.unlock();
    
    return true;
    
}

InRead* InReads::traverseInRead(Log* threadLog, Sequence* sequence, unsigned int seqPos) { // traverse a single read

    unsigned long long int A = 0, C = 0, G = 0, T = 0, N = 0, lowerCount = 0;
    unsigned long long int sumQuality =0;
    double avgQuality=0;
    
    for (char &base : *sequence->sequence) {
        
        if (islower(base)) {
            
            lowerCount++;
            
        }
                
        switch (base) {
            case 'A':
            case 'a':{
                
                ++A;
                break;
                
            }
            case 'C':
            case 'c':{
                
                ++C;
                break;
                
            }
            case 'G':
            case 'g': {
                
                ++G;
                break;
                
            }
            case 'T':
            case 't': {
                
                ++T;
                break;
                
            }

            case 'N':
            case 'n':
            case 'X':
            case 'x': {

                ++N;
                break;
            }
                
            default: {
                break;
            }
                
        }
    }


    if (sequence->sequenceQuality != NULL){
    for (char &quality : *sequence -> sequenceQuality) { 

        sumQuality += int(quality) - 33;

        }
    
    avgQuality = sumQuality/(sequence->sequenceQuality->size());
    }

        // operations on the segment

    InRead* inRead = new InRead;
    
    inRead->set(threadLog, 0, 0, sequence->header, &sequence->comment, sequence->sequence, &A, &C, &G, &T, &lowerCount, seqPos, sequence->sequenceQuality, &avgQuality, NULL, &N);

    return inRead;
    
}

unsigned long long int InReads::getTotReadLen() {
    
    // unsigned long long int totReadLen;
    
    // for (InRead* read : inReads) {
        
    //     totReadLen = read->getA() + read->getC() + read->getG() + read->getT();
    //     std::cout << totReadLen << "\n";
        
    // }
    
    return totA + totC + totG + totT + totN;
    
}

int InReads::computeGCcontent() { 

    unsigned long long int totReadLen = totA + totC + totG + totT; 
    double GCcontent = (double) (totG+totC)/totReadLen * 100;
    
    return GCcontent;
}

double InReads::computeAvgReadLen() {
    
    return (double) getTotReadLen()/inReads.size();
    
}

unsigned long long int InReads::getReadN50() {
    
    return readNstars[4];
    
}

void InReads::evalNstars() {

    computeNstars(readLens, readNstars, readLstars); 
    
}

int InReads::getSmallestRead() {

    return readLens.back();

}

int InReads::getLargestRead() {

    return readLens.front();

}

void InReads::getQualities(){

    for (const auto &item : qualities){
        std::cout << item << "\n";
    }
    std::cout << std::endl;


}


double InReads::getAvgQualities(){

unsigned long long int sumQualities = 0;
unsigned long long int avgQualitiesSize=avgQualities.size();
double avgQuality = 0;

for (unsigned long long int i=0; i < avgQualitiesSize; i++) {

    sumQualities += avgQualities[i]*readLens[i];  //sum the qualities normalized by their read length

}

avgQuality = (sumQualities/getTotReadLen());
return avgQuality;

}


void InReads::report(unsigned long long int gSize) {

    if (inReads.size() > 0) {
        
        if (!tabular_flag) {
        
            std::cout<<output("+++Read summary+++")<<"\n";
        
        }
        
        std::cout<<output("# reads")<<inReads.size()<<"\n";
        std::cout<<output("Total read length")<<getTotReadLen()<<"\n";
        std::cout<<output("Average read length") << gfa_round(computeAvgReadLen()) << "\n";
        evalNstars(); // read N* statistics
        std::cout<<output("Read N50")<<getReadN50()<<"\n";
        std::cout<<output("Smallest read length")<<getSmallestRead()<<"\n";
        std::cout<<output("Largest read length")<<getLargestRead()<<"\n";
        std::cout<<output("Coverage")<<gfa_round((double)getTotReadLen()/gSize)<<"\n";
        std::cout<<output("GC content %")<<computeGCcontent()<<"\n";
        std::cout<<output("Base composition (A:C:T:G)")<<totA<<":"<<totC<<":"<<totT<<":"<<totG<<"\n";
        std::cout<<output("Average read quality")<<getAvgQualities()<<"\n";
        
    }
    
}

// bool compareReadLengths(InRead* read1, InRead* read2) {
//     return (read1->getReadLen() < read2->getReadLen());

// }

void InReads::printReadLengths(char sizeOutType) {

    if (sizeOutType == 's' || sizeOutType == 'h' || sizeOutType == 'c') {
        sort(readLens.begin(), readLens.end());

    }

    if (sizeOutType == 'h') {

        int count = 1; 
        for (unsigned long long int i = 0; i < readLens.size(); i++) {
            if (readLens[i] == readLens[i+1]) {
                count += 1;
            }
            else if (readLens[i] != readLens[i+1]) {
                std::cout << readLens[i] << "," << count << "\n";
                count = 1;
            }
        }
    }

    if (sizeOutType == 'c') {

        int count = 1; 
        unsigned long long int sizexCount;
        std::vector<unsigned int> counts;
        std::vector<unsigned long long int> sizexCounts;
        unsigned long long int sizexCountSum = 0;
        std::vector<unsigned long long int> uniqReadLens;

        for (unsigned long long int i = 0; i < readLens.size(); i++) {
            if (readLens[i] == readLens[i+1]) {
                count += 1; 
            }
            else if (readLens[i] != readLens[i+1]) {
                sizexCount = (readLens[i] * count);
                counts.push_back(count);
                sizexCounts.push_back(sizexCount);
                uniqReadLens.push_back(readLens[i]);
                sizexCountSum += sizexCount;

                count = 1;
            }
        }

        unsigned long long int sizexCountSums = 0; 
        for (unsigned long long int i = 0; i < sizexCounts.size(); i++) {
            std::cout << uniqReadLens[i] << "," << counts[i] << "," <<sizexCounts[i] << "," << sizexCountSum - sizexCountSums << "\n"; 
            sizexCountSums += sizexCounts[i];
        }


    }

    else {
        for (auto i: readLens) {
            std::cout << i << "\n";
        }
    }

}

void InReads::printQualities(char qualityOut) {

    if (qualityOut == 'c'){ 
        for (unsigned long long int i = 0; i < (avgQualities.size()); i++) {
            std::cout << avgQualities[i] << "\n"; 
        }
    }
    else if (qualityOut == 'l') { // l prints read lengths and qualities 
        for (unsigned long long int i = 0; i < (avgQualities.size()); i++) {
            std::cout << readLens[i] << "," << avgQualities[i] << "\n";
        }
    }

}

InRead::~InRead()
{
    delete inRead;
    delete inSequenceQuality;
}

void InRead::set(Log* threadLog, unsigned int uId, unsigned int iId, std::string seqHeader, std::string* seqComment, std::string* sequence, unsigned long long int* A, unsigned long long int* C, unsigned long long int* G, unsigned long long int* T, unsigned long long int* lowerCount, unsigned int seqPos, std::string* sequenceQuality, double* avgQuality, std::vector<Tag>* inSequenceTags, unsigned long long int* N) {
    
    threadLog->add("Processing read: " + seqHeader + " (uId: " + std::to_string(uId) + ", iId: " + std::to_string(iId) + ")");
    
    unsigned long long int seqSize = 0;
    
    this->setiId(iId); // set temporary sId internal to scaffold
    
    this->setuId(uId); // set absolute id
    
    this->setReadPos(seqPos); // set original order
    
    this->setReadHeader(&seqHeader);
    
    if (*seqComment != "") {
        
        this->setReadComment(*seqComment);
        
    }
    
    if (inSequenceTags != NULL) {
        
        this->setReadTags(inSequenceTags);
        
    }
    
    if (*sequence != "*") {
        
        this->setInRead(sequence);
        
        threadLog->add("Segment sequence set");
        
        if (sequenceQuality != NULL) {
            
            this->setInReadQuality(sequenceQuality);
            
            threadLog->add("Segment sequence quality set");

            this->setAvgQuality(avgQuality);
            
        }
        
        this->setACGT(A, C, G, T, N);
        
        threadLog->add("Increased ACGT counts");
        
        this->setLowerCount(lowerCount);

        threadLog->add("Increased total count of lower bases");
        
        seqSize = *A + *C + *G + *T;
        
    }else{
        
        seqSize = *lowerCount;
        
        this->setLowerCount(&seqSize);
        
        threadLog->add("No seq input. Length (" + std::to_string(seqSize) + ") recorded in lower count");
        
    }
    
}
    
void InRead::setReadHeader(std::string* h) {
    seqHeader = *h;
}

void InRead::setReadComment(std::string c) {
    seqComment = c;
}

void InRead::setInRead(std::string* s) {
    inRead = s;
}

void InRead::setInReadQuality(std::string* q) {
    inSequenceQuality = q;
}

void InRead::setReadTags(std::vector<Tag>* tag) {
    tags = *tag;
}


void InRead::setuId(unsigned int i) { // absolute id
    uId = i;
}

void InRead::setiId(unsigned int i) { // temporary id, internal to scaffold
    iId = i;
}

void InRead::setReadPos(unsigned int i) { // temporary id, internal to scaffold
    seqPos = i;
}

std::string InRead::getReadHeader() {
    return seqHeader;
}

std::string InRead::getReadComment() {
    return seqComment;
}

std::vector<Tag> InRead::getTags() {
    return tags;
}

std::string InRead::getInRead(unsigned int start, unsigned int end) {
    
    if (inRead == NULL) {
        
        return "*";
        
    }else{
    
        return start != 0 || end != 0 ? inRead->substr(start-1, end-start+1) : *inRead;
        
    }
    
}

std::string InRead::getInReadQuality(unsigned int start, unsigned int end) {
    
    if (inSequenceQuality != NULL) {
    
        return start != 0 || end != 0 ? inSequenceQuality->substr(start-1, end-start+1) : *inSequenceQuality;
        
    }else{
        
        return "";
        
    }
    
}

// std::vector<int> InRead::getQualitiesInt() {

//     return QualitiesInt;

// }

unsigned int InRead::getReadPos() {
    
    return seqPos;

}

unsigned long long int InRead::getReadLen(unsigned long long int start, unsigned long long int end) {
    
    if (inRead == NULL) {
        
        return lowerCount;
        
    }else{
    
        return start != 0 || end != 0 ? end-start+1 : A + C + G + T + N; // need to sum long long int to prevent size() overflow
        
    }
    
}

unsigned int InRead::getuId() { // absolute id
    
    return uId;
}

unsigned int InRead::getiId() { // temporary id, internal to scaffold
    
    return iId;
}

void InRead::setAvgQuality(double* InReadAvgQuality) {
    
    avgQuality = *InReadAvgQuality;

}

// void InRead::setQualitiesInt(std::vector<int>* qualInt) {

//     QualitiesInt = *qualInt;

// }

void InRead::setACGT(unsigned long long int* a, unsigned long long int* c, unsigned long long int* g, unsigned long long int* t, unsigned long long int* n) {
    
    A = *a;
    C = *c;
    G = *g;
    T = *t;
    N = *n;
    
}

void InRead::setLowerCount(unsigned long long int* C) {
    
    lowerCount = *C;
    
}

unsigned long long int InRead::getA() {
    
    return A;
}

unsigned long long int InRead::getC() {
    
    return C;
}

unsigned long long int InRead::getG() {
    
    return G;
}

unsigned long long int InRead::getT() {
    
    return T;
}

unsigned long long int InRead::getN() {

    return N;
}



unsigned int InRead::getLowerCount(unsigned long long int start, unsigned long long int end) {
    
    if (start == 0 || end == 0) {
        
        return lowerCount;
        
    }else{
        
        unsigned long long int lowerCountSubset = 0;
        
        for (char base : *inRead) { // need to fix this loop
            
            if (islower(base)) {
                
                ++lowerCountSubset;
                
            }
            
        }
        
        return lowerCountSubset;
        
    }

}

double InRead::computeGCcontent() {
    
    double GCcontent = (double) (G + C) / (G + C + A + T + N) * 100;
    
    return GCcontent;
}

double InRead::getAvgQuality() {

    return avgQuality;

}

bool InRead::trimRead(unsigned int start, unsigned int end) {
    
    for(char& base : inRead->substr(start, end-start)) {
        
        switch (base) {
            case 'A':
            case 'a':{
                
                A--;
                break;
                
            }
            case 'C':
            case 'c':{
                
                C--;
                break;
                
            }
            case 'G':
            case 'g': {
                
                G--;
                break;
                
            }
            case 'T':
            case 't': {
                
                T--;
                break;
                
            }
                
        }
        
    }
    
    inRead->erase(start, end-start);
    
    if (inSequenceQuality->size()>0) {
    
        inSequenceQuality->erase(start, end-start);
    
    }
    
    return true;
}

bool InRead::rvcpRead() {

    *inRead = revCom(*inRead);

    return true;
    
}

bool InRead::invertRead() {

    *inRead = rev(*inRead);
    
    if (inSequenceQuality != NULL) {
    
        *inSequenceQuality = rev(*inSequenceQuality);
    
    }
        
    return true;
    
}
