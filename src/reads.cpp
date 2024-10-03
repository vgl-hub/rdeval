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
#include "uid-generator.h"
#include "gfa.h"
#include "functions.h" // global functions
#include "stream-obj.h"

#include "zlib.h"
#include "zstream/zstream_common.hpp"
#include "zstream/ozstream.hpp"
#include "zstream/ozstream_impl.hpp"
#include "output.h"

#include "reads.h"

void InRead::set(Log* threadLog, uint32_t uId, uint32_t iId, std::string seqHeader, std::string* seqComment, std::string* sequence, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, uint32_t seqPos, std::string* sequenceQuality, double* avgQuality, std::vector<Tag>* inSequenceTags, uint64_t* N) {
    
    threadLog->add("Processing read: " + seqHeader + " (uId: " + std::to_string(uId) + ", iId: " + std::to_string(iId) + ")");
    uint64_t seqSize = 0;
    this->setiId(iId); // set temporary sId internal to scaffold
    this->setuId(uId); // set absolute id
    this->setSeqPos(seqPos); // set original order
    this->setSeqHeader(seqHeader);
    if (*seqComment != "")
        this->setSeqComment(*seqComment);
    
    if (inSequenceTags != NULL)
        this->setSeqTags(inSequenceTags);
    
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

void InReads::load() {
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    std::size_t numFiles = userInput.inFiles.size();
    lg.verbose("Processing " + std::to_string(numFiles) + " files");
    
    for (uint32_t i = 0; i < numFiles; i++) {
        
        std::string ext = getFileExt(userInput.file('r', i));
        
        if (ext == "bam") {
            
        }else if (ext == "fasta" || ext == "fastq" || ext == "fasta.gz" || ext == "fastq.gz") {
            
            StreamObj streamObj;
            stream = streamObj.openStream(userInput, 'r', i);
            Sequences* readBatch = new Sequences;
            
            if (stream) {
                
                switch (stream->peek()) {
                        
                    case '>': {
                        
                        stream->get();
                        
                        while (getline(*stream, newLine)) {
                            
                            h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                            c = strtok(NULL,""); //read comment
                            
                            seqHeader = h;
                            
                            if (c != NULL)
                                seqComment = std::string(c);
                            
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
                            
                            if (c != NULL)
                                seqComment = std::string(c);
                            
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
        }else{
            fprintf(stderr, "cannot recognize input (must be: fasta, fastq, bam, cram).\n");
            exit(EXIT_FAILURE);
        }
    }
}

void InReads::appendReads(Sequences* readBatch) { // read a collection of reads
    
    threadPool.queueJob([=]{ return traverseInReads(readBatch); });
    writeToStream();
}

bool InReads::traverseInReads(Sequences* readBatch) { // traverse the read

    Log threadLog;
    threadLog.setId(readBatch->batchN);
    std::vector<InRead*> inReadsBatch;
    uint32_t readN = 0;
    std::vector<uint64_t> readLensBatch;
    InRead* read;

    uint64_t batchA = 0, batchT=0, batchC=0, batchG=0, batchN =0;
    // std::vector<long double> batchListA, batchListC, batchListT, batchListG, batchListN;
    std::vector<double> batchAvgQualities;
    uint64_t filterInt = 0;

    if (userInput.filter != "none")
        filterInt = stoi(userInput.filter.substr(1));
    
    for (Sequence* sequence : readBatch->sequences) {

        if (userInput.filter != "none" &&
            (
                ((userInput.filter[0] == '>') && (sequence->sequence->size() <= filterInt)) ||
                ((userInput.filter[0] == '<') && (sequence->sequence->size() >= filterInt)) ||
                ((userInput.filter[0] == '=') && (sequence->sequence->size() != filterInt))
            )
        ){
            
            lg.verbose("Sequence length (" + std::to_string(sequence->sequence->size()) + ") shorter than filter length. Filtering out (" + sequence->header + ").");
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
        
        inReadsBatch.push_back(read);
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    readBatches.emplace_back(inReadsBatch,readBatch->batchN);
    delete readBatch;
    readLens.insert(std::end(readLens), std::begin(readLensBatch), std::end(readLensBatch));
    avgQualities.insert(std::end(avgQualities), std::begin(batchAvgQualities), std::end(batchAvgQualities));
    totReads += inReadsBatch.size();

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
        for (char &quality : *sequence -> sequenceQuality)
            sumQuality += int(quality) - 33;
        avgQuality = sumQuality/(sequence->sequenceQuality->size());
    }

    // operations on the segment
    InRead* inRead = new InRead;
    inRead->set(threadLog, 0, 0, sequence->header, &sequence->comment, sequence->sequence, &A, &C, &G, &T, &lowerCount, seqPos, sequence->sequenceQuality, &avgQuality, NULL, &N);
    sequence->sequence = NULL;
    sequence->sequenceQuality = NULL;
    return inRead;
}

uint64_t InReads::getTotReadLen() {
    return totA + totC + totG + totT + totN;
}

double InReads::computeGCcontent() {

    uint64_t totReadLen = totA + totC + totG + totT;
    double GCcontent = (double) (totG+totC)/totReadLen * 100;
    
    return GCcontent;
}

double InReads::computeAvgReadLen() {
    return (double) getTotReadLen()/totReads;
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

void InReads::report() {

    if (totReads > 0) {
        
        if (!tabular_flag)
            std::cout<<output("+++Read summary+++")<<"\n";
        std::cout<<output("# reads")<<totReads<<"\n";
        std::cout<<output("Total read length")<<getTotReadLen()<<"\n";
        std::cout<<output("Average read length") << gfa_round(computeAvgReadLen()) << "\n";
        evalNstars(); // read N* statistics
        std::cout<<output("Read N50")<<getReadN50()<<"\n";
        std::cout<<output("Smallest read length")<<getSmallestRead()<<"\n";
        std::cout<<output("Largest read length")<<getLargestRead()<<"\n";
        std::cout<<output("Coverage")<<gfa_round((double)getTotReadLen()/userInput.gSize)<<"\n";
        std::cout<<output("GC content %")<<computeGCcontent()<<"\n";
        std::cout<<output("Base composition (A:C:T:G)")<<totA<<":"<<totC<<":"<<totT<<":"<<totG<<"\n";
        std::cout<<output("Average read quality")<<getAvgQualities()<<"\n";
    }
}

void InReads::printReadLengths() {

    if (userInput.sizeOutType == 's' || userInput.sizeOutType == 'h' || userInput.sizeOutType == 'c') {
        sort(readLens.begin(), readLens.end());

    }

    if (userInput.sizeOutType == 'h') {

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

    if (userInput.sizeOutType == 'c') {

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

void InReads::printQualities() {
    
    char qualityOut = userInput.qualityOut;

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
void InReads::printContent() {
    
    char content = userInput.content;

    for (std::pair<std::vector<InRead*>,uint32_t> inReads : readBatches){
        
        for (InRead* read : inReads.first){
            
            long double readLen=read->getA()+read->getT()+read->getC()+read->getG()+read->getN();
            
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
}

void InReads::writeToStream() {
    
    if (userInput.outFiles.size()) {
        
        std::string ext = getFileExt(outputStream.file); // variable to handle output path and extension
        
        const static phmap::flat_hash_map<std::string,int> string_to_case{
            {"fasta",1},
            {"fa",1},
            {"fasta.gz",1},
            {"fa.gz",1},
            {"fastq",2},
            {"fq",2},
            {"fastq.gz",2},
            {"fq.gz",2},
            
        };
        
        std::vector<std::pair<std::vector<InRead*>,uint32_t>> readBatchesCpy;
        {
            std::unique_lock<std::mutex> lck(mtx);
            readBatchesCpy = {readBatches.begin() + batchCounter-1, readBatches.end()};
        }
        
        switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
                
            case 1:  {// fasta[.gz]
                
                for (std::pair<std::vector<InRead*>,uint32_t> inReads : readBatchesCpy) {
                    
                    if (inReads.second > batchCounter)
                        continue;
                    
                    lg.verbose("Writing read batch " + std::to_string(inReads.second) + " to file (" + std::to_string(inReads.first.size())  + ")");
                    
                    for (InRead* read : inReads.first){
                        
                        *outputStream.stream << '>' << read->seqHeader << '\n' << *read->inSequence << '\n';
                        delete read;
                    }
                    ++batchCounter;
                }
                break;
            }
            
            case 2:  {// fastq[.gz]
                
                for (std::pair<std::vector<InRead*>,uint32_t> inReads : readBatchesCpy) {
                    
                    if (inReads.second > batchCounter)
                        continue;
                    
                    lg.verbose("Writing read batch " + std::to_string(inReads.second) + " to file (" + std::to_string(inReads.first.size())  + ")");
                        
                        for (InRead* read : inReads.first){
                            
                            *outputStream.stream << '@' << read->seqHeader << '\n' << *read->inSequence << "\n+\n" << (read->inSequenceQuality != NULL ? *read->inSequenceQuality : std::string('!', read->inSequence->size())) << '\n';
                            delete read;
                        }
                    ++batchCounter;
                }
                break;
            }
        }
    }
}
