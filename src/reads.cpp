#include <stdlib.h>
#include <string>
#include <vector>
#include <mutex>
#include <string.h>

#include <iostream>
#include <fstream>

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

#include "log.h"
#include "global.h"

#include "bed.h"
#include "struct.h"
#include "gfa-lines.h"
#include "stream-obj.h"

#include "functions.h" // global functions

#include "reads.h"

InReads::~InReads()
{

    for (InSegment* p : inReads)
        delete p;

}

void InReads::load(UserInput userInput) {

    unsigned int batchSize = 100;
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    std::shared_ptr<std::istream> stream;
    
    unsigned int numFiles = userInput.iReadFileArg.size();
    
    lg.verbose("Processing " + std::to_string(numFiles) + " files");
    
    for (unsigned int i = 0; i < numFiles; i++) {
        
        StreamObj streamObj;

        stream = streamObj.openStream(userInput, 'r', &i);

        Sequences* readBatch = new Sequences;

        if (stream) {

            switch (stream->peek()) {

                case '>': {

                    stream->get();

                    while (!stream->eof()) {

                        getline(*stream, newLine);

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

                            appendReads(readBatch);

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

                        }else{

                            seqComment = "";

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

                            appendReads(readBatch);

                            readBatch = new Sequences;

                        }

                        lg.verbose("Individual fastq sequence read: " + seqHeader);

                    }

                    break;

                }

            }
            
            readBatch->batchN = seqPos/batchSize + 1;
                
            lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));

            appendReads(readBatch);

        }
        
    }

}

void InReads::appendReads(Sequences* readBatch) { // read a collection of reads
    
    threadPool.queueJob([=]{ return traverseInReads(readBatch); });
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
    lck.unlock();
    
}

bool InReads::traverseInReads(Sequences* readBatch) { // traverse the read

    Log threadLog;
    
    threadLog.setId(readBatch->batchN);
    
    std::vector<InSegment*> inReadsBatch;
    
    unsigned int readN = 0; 
    std::vector<unsigned long long int> readLensBatch;
    InSegment* read;
    unsigned long long int batchA = 0, batchT=0, batchC=0, batchG=0;

    for (Sequence* sequence : readBatch->sequences) {
        
        read = traverseInRead(&threadLog, sequence, readBatch->batchN+readN++);
        readLensBatch.push_back(read->getSegmentLen());
        batchA += read->getA();
        batchT += read->getT();
        batchC += read->getC();
        batchG += read->getG();

        inReadsBatch.push_back(read);
        
    }
    
    delete readBatch;
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    inReads.insert(std::end(inReads), std::begin(inReadsBatch), std::end(inReadsBatch));
    readLens.insert(std::end(readLens), std::begin(readLensBatch), std::end(readLensBatch));
    totA+=batchA;
    totT+=batchT;
    totC+=batchC;
    totG+=batchG;

    logs.push_back(threadLog);
    
    lck.unlock();
    
    return true;
    
}

InSegment* InReads::traverseInRead(Log* threadLog, Sequence* sequence, unsigned int seqPos) { // traverse a single read
    
    unsigned long long int A = 0, C = 0, G = 0, T = 0, N = 0, lowerCount = 0;
    
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
    
    // operations on the segment
    InSegment* inSegment = new InSegment;
    
    inSegment->set(threadLog, 0, 0, sequence->header, &sequence->comment, sequence->sequence, &A, &C, &G, &T, &lowerCount, seqPos, sequence->sequenceQuality, NULL, &N);
    
    return inSegment;
    
}

unsigned long long int InReads::getTotReadLen() {
    
    return totA + totC + totG + totT;
    
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

    computeNstars(readLens, readNstars, readLstars); // WOULD BE HELPFUL TO READLENS OUTSIDE OF THIS FUNCTION
    
}

int InReads::getSmallestRead() {

    return readLens.back();

}

int InReads::getLargestRead() {

    return readLens.front();

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


bool compareReadLengths(InRead* read1, InRead* read2) {
    return (read1->getReadLen() < read2->getReadLen());

}

void InReads::printReadLengths(char sizeOutType) {

    // if (sizeOutType == 's') {
    //     sort(inReads.begin(), inReads.end(), compareReadLengths);
    // }

    // for (InRead* read : inReads) {
        
    //     std::cout << (read->getA() + read->getC() + read->getG() + read->getT()) << "\n";
        
    // }
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

        
    }
    
}
