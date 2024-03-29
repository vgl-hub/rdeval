#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fstream>
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

void InRead::set(Log* threadLog, uint32_t uId, uint32_t iId, std::string seqHeader, std::string* seqComment, std::string* sequence, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, uint32_t seqPos, std::string* sequenceQuality, double* avgQuality, std::vector<Tag>* inSequenceTags, uint64_t* N) {
    
    threadLog->add("Processing read: " + seqHeader + " (uId: " + std::to_string(uId) + ", iId: " + std::to_string(iId) + ")");
    
    uint64_t seqSize = 0;
    
    this->setiId(iId); // set temporary sId internal to scaffold
    
    this->setuId(uId); // set absolute id
    
    this->setSeqPos(seqPos); // set original order
    
    this->setSeqHeader(&seqHeader);
    
    if (*seqComment != "") {
        
        this->setSeqComment(*seqComment);
        
    }
    
    if (inSequenceTags != NULL) {
        
        this->setSeqTags(inSequenceTags);
        
    }
    
    if (*sequence != "*") {
        
        this->setInSequence(sequence);
        
        threadLog->add("Segment sequence set");
        
        if (sequenceQuality != NULL) {
            
            this->setInSequenceQuality(sequenceQuality);
            
            threadLog->add("Segment sequence quality set");

            this->avgQuality = *avgQuality;
            
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

void InReads::load(UserInputRdeval* userInput) {
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    uint32_t numFiles = userInput->inReads.size();
    
    lg.verbose("Processing " + std::to_string(numFiles) + " files");
    
    for (uint32_t i = 0; i < numFiles; i++) {
        
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
    
    std::unique_lock<std::mutex> lck (mtx);
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
}

bool InReads::traverseInReads(Sequences* readBatch, UserInputRdeval* userInput) { // traverse the read

    Log threadLog;
    
    threadLog.setId(readBatch->batchN);
    
    std::vector<InRead*> inReadsBatch;
    
    uint32_t readN = 0;
    std::vector<uint64_t> readLensBatch;
    InRead* read;
    uint64_t batchA = 0, batchT=0, batchC=0, batchG=0, batchN =0;
    // std::vector<long double> batchListA, batchListC, batchListT, batchListG, batchListN;
    std::vector<double> batchAvgQualities;
    uint32_t filterInt = 0;
    if (!(userInput->filter == "none")) {
        filterInt = stoi(userInput->filter.substr(1));
    }


    // perhaps I need to hard code a small function because the filter is currently in the for-loop and basically requires re-parsing the filter operand each time

    for (Sequence* sequence : readBatch->sequences) {

        if ((userInput->filter[0] == '>') && (sequence->sequence->size() <= filterInt)){
            // std::cout << "Sequence shorter than filter length." << "\n"; // less compute time to assign to a variable or call like this? // would like to output sequence name here 
            continue;
        }
        else if ((userInput->filter[0] == '<') && (sequence->sequence->size() >= filterInt)) {
            // std::cout << "Sequence longer than filter length." << "\n"; 
            continue;
        }
        else if ((userInput->filter[0] == '=') && (sequence->sequence->size() != filterInt)) {
            // std::cout << "Sequence does not equal filter length." << "\n";
            continue;
        }
        
        read = traverseInRead(&threadLog, sequence, readBatch->batchN+readN++);

        readLensBatch.push_back(read->inSequence->size());
        batchA += read->getA();
        batchT += read->getT();
        batchC += read->getC();
        batchG += read->getG();
        batchN += read->getN();

        if (read->inSequenceQuality != NULL) 
            batchAvgQualities.push_back(read->avgQuality);

        // // check for flag such that if there isn't a flag it doesn't create these 
        // batchListA.push_back(read -> getA()/float(sequence->sequence->size()));
        // //batchListA.push_back(read->getA());
        // batchListT.push_back(read -> getT()/float(sequence->sequence->size()));
        // batchListC.push_back(read -> getC()/float(sequence->sequence->size()));
        // batchListG.push_back(read -> getG()/float(sequence->sequence->size()));
        // batchListN.push_back(read -> getN()/float(sequence->sequence->size()));
        
        inReadsBatch.push_back(read);
        
    }
    
    delete readBatch;
    
    std::unique_lock<std::mutex> lck(mtx);
    
    inReads.insert(std::end(inReads), std::begin(inReadsBatch), std::end(inReadsBatch));
    readLens.insert(std::end(readLens), std::begin(readLensBatch), std::end(readLensBatch));
    avgQualities.insert(std::end(avgQualities), std::begin(batchAvgQualities), std::end(batchAvgQualities));
    // listA.insert(std::end(listA),std::begin(batchListA),std::end(batchListA));

    totA+=batchA;
    totT+=batchT;
    totC+=batchC;
    totG+=batchG;
    totN+=batchN;

    logs.push_back(threadLog);
    
    return true;
    
}

InRead* InReads::traverseInRead(Log* threadLog, Sequence* sequence, uint32_t seqPos) { // traverse a single read

    uint64_t A = 0, C = 0, G = 0, T = 0, N = 0, lowerCount = 0;
    uint64_t sumQuality =0;
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

uint64_t InReads::getTotReadLen() {
    
    // uint64_t totReadLen;
    
    // for (InRead* read : inReads) {
        
    //     totReadLen = read->getA() + read->getC() + read->getG() + read->getT();
    //     std::cout << totReadLen << "\n";
        
    // }
    
    return totA + totC + totG + totT + totN;
    
}

double InReads::computeGCcontent() {

    uint64_t totReadLen = totA + totC + totG + totT;
    double GCcontent = (double) (totG+totC)/totReadLen * 100;
    
    return GCcontent;
}

double InReads::computeAvgReadLen() {
    
    return (double) getTotReadLen()/inReads.size();
    
}

uint64_t InReads::getReadN50() {
    
    return readNstars[4];
    
}

void InReads::evalNstars() {

    computeNstars(readLens, readNstars, readLstars); 
    
}

uint64_t InReads::getSmallestRead() {

    return readLens.back();

}

uint64_t InReads::getLargestRead() {

    return readLens.front();

}

void InReads::getQualities(){

    for (const auto &item : qualities)
        std::cout << item << "\n";
    
    std::cout << std::endl;


}


double InReads::getAvgQualities(){

    uint64_t sumQualities = 0;
    uint64_t avgQualitiesSize=avgQualities.size();
    double avgQuality = 0;

    for (uint64_t i=0; i < avgQualitiesSize; i++)
        sumQualities += avgQualities[i]*readLens[i];  //sum the qualities normalized by their read length

    avgQuality = (sumQualities/getTotReadLen());
    return avgQuality;

}


void InReads::report(uint64_t gSize) {

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
        for (uint64_t i = 0; i < readLens.size(); i++) {
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
        uint64_t sizexCount;
        std::vector<unsigned  int> counts;
        std::vector<uint64_t> sizexCounts;
        uint64_t sizexCountSum = 0;
        std::vector<uint64_t> uniqReadLens;

        for (uint64_t i = 0; i < readLens.size(); i++) {
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

        uint64_t sizexCountSums = 0;
        for (uint64_t i = 0; i < sizexCounts.size(); i++) {
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
        for (uint64_t i = 0; i < (avgQualities.size()); i++) {
            std::cout << avgQualities[i] << "\n"; 
        }
    }
    else if (qualityOut == 'l') { // l prints read lengths and qualities 
        for (uint64_t i = 0; i < (avgQualities.size()); i++) {
            std::cout << readLens[i] << "," << avgQualities[i] << "\n";
        }
    }

}
void InReads::printContent(char content) {
    // if (content == 'a'){ //imaging a function in which you could perhaps print different combinations of bases, or just N's or in different formats (i.e./ percent of reads v. normalized v. not etc.)
    //     for (uint64_t i = 0; i < (listA.size()); i++) {
    //         std::cout << listA[i] << "\n";
    //     }

    // }

    for (InRead* read : inReads){

        long double readLen=read->getA()+read->getT()+read->getC()+read->getG()+read->getN();

        // if (content == 'a' && outflag == ){
        //     std::cout << read->getA()/readLen << "," << read->getT()/readLen << "," << read->getC()/readLen << "," << read->getG()/readLen <<"," << read->getN()/readLen << "\n";
        // }

        if (content == 'g') {
            std::cout << (read->getC()+read->getG())/readLen << "\n";
        }

        if (content == 't') {
            std::cout << (read->getA()+read->getT())/readLen << "\n";
        }

        if (content == 'n') {
            std::cout << read->getN()/readLen << "\n";
        }

    }
}
