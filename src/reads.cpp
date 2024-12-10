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
#include <iomanip>

#include "zlib.h"
#include "zstream/zstream_common.hpp"
#include "zstream/ozstream.hpp"
#include "zstream/ozstream_impl.hpp"
#include "output.h"
#include "len-vector.h"

#include "reads.h"

void InRead::set(Log* threadLog, uint32_t uId, uint32_t iId, std::string seqHeader, std::string* seqComment, std::string* sequence, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, uint32_t seqPos, std::string* sequenceQuality, double avgQuality, std::vector<Tag>* inSequenceTags, uint64_t* N) {
    
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
            this->avgQuality = avgQuality;
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
    uint32_t batchN = 0;
    uint64_t processedLength = 0;
    lg.verbose("Processing " + std::to_string(numFiles) + " files");
    
    const static phmap::flat_hash_map<std::string,int> string_to_case{
        {"fasta",1},
        {"fa",1},
        {"fasta.gz",1},
        {"fa.gz",1},
        {"fastq",1},
        {"fq",1},
        {"fastq.gz",1},
        {"fq.gz",1},
        {"bam",2},
        {"cram",3},
        {"rd",4}
    };
    
    for (uint32_t i = 0; i < numFiles; i++) {
        
        std::string ext = getFileExt(userInput.file('r', i));
        
        switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
                
            case 1: { // fa*[.gz]
                
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
                                processedLength += inSequence->size();
                                
                                if (processedLength > batchSize) {
                                    readBatch->batchN = ++batchN;
                                    lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
                                    appendReads(readBatch);
                                    readBatch = new Sequences;
                                    processedLength = 0;
                                }
                                //lg.verbose("Individual fasta sequence read: " + seqHeader);
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
                                processedLength += inSequence->size();
                                
                                if (processedLength > batchSize) {
                                    readBatch->batchN = ++batchN;
                                    lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
                                    appendReads(readBatch);
                                    readBatch = new Sequences;
                                    processedLength = 0;
                                }
                                //lg.verbose("Individual fastq sequence read: " + seqHeader);
                            }
                            break;
                        }
                    }
                    readBatch->batchN = ++batchN; // process residual reads
                    lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
                    appendReads(readBatch);
                }
                break;
            }
            case 2: { // bam
                
                fprintf(stderr, "bam currently not supported.\n");
                exit(EXIT_FAILURE);
                
                StreamObj streamObj;
                stream = streamObj.openStream(userInput, 'r', i);
                //Sequences* readBatch = new Sequences;
                
                while (getline(*stream, newLine)) {
                    
                    std::cout<<newLine<<std::endl;
                    
                }
                break;
            }
            case 3: { // cram
                fprintf(stderr, "cram currently not supported.\n");
                exit(EXIT_FAILURE);
                break;
            }
            case 4: { // rd
                    readTableCompressedBinary(userInput.inFiles[i]);
                break;
            }
            default: { // fasta[.gz]
                fprintf(stderr, "cannot recognize input (must be: fasta, fastq, bam, cram).\n");
                exit(EXIT_FAILURE);
            }
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
    LenVector<float> readLensBatch;
    InRead* read;

    uint64_t batchA = 0, batchT=0, batchC=0, batchG=0, batchN =0;
    std::vector<float> batchAvgQualities;
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
        std::pair<uint64_t, float> lenQual(read->inSequence->size(), read->avgQuality);
        readLensBatch.push_back(lenQual);
        batchA += read->getA();
        batchT += read->getT();
        batchC += read->getC();
        batchG += read->getG();
        batchN += read->getN();

        if (read->inSequenceQuality != NULL)
            batchAvgQualities.push_back(read->avgQuality);
        
        if (streamOutput)
            inReadsBatch.push_back(read);
        else
            delete read;
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    readBatches.emplace_back(inReadsBatch,readBatch->batchN);
    delete readBatch;
    readLens.insert(readLensBatch);
    avgQualities.insert(std::end(avgQualities), std::begin(batchAvgQualities), std::end(batchAvgQualities));

    totA+=batchA;
    totT+=batchT;
    totC+=batchC;
    totG+=batchG;
    totN+=batchN;
    totReads += readN;

    logs.push_back(threadLog);
    return true;
}

InRead* InReads::traverseInRead(Log* threadLog, Sequence* sequence, uint32_t seqPos) { // traverse a single read

    uint64_t A = 0, C = 0, G = 0, T = 0, N = 0, lowerCount = 0;
    uint64_t sumQuality = 0;
    float avgQuality = 0;
    
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
        avgQuality = (float) sumQuality/(sequence->sequenceQuality->size());
    }

    // operations on the segment
    InRead* inRead = new InRead;
    inRead->set(threadLog, 0, 0, sequence->header, &sequence->comment, sequence->sequence, &A, &C, &G, &T, &lowerCount, seqPos, sequence->sequenceQuality, avgQuality, NULL, &N);
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

void InReads::evalNstars() { // not very efficient
    std::vector<std::pair<uint64_t,float>> allReadLens = readLens.all();
    std::vector<uint64_t> tmpVector(allReadLens.size());
    for (uint64_t i = 0; i < allReadLens.size(); ++i)
        tmpVector[i] = allReadLens[i].first;
    computeNstars(tmpVector, readNstars, readLstars);
}

uint64_t InReads::getSmallestRead() {
    return readLens.front();
}

uint64_t InReads::getLargestRead() {
    return readLens.back();
}

double InReads::getAvgQuality(){

    uint64_t sumQualities = 0, avgQualitiesSize=avgQualities.size();

    for (uint64_t i = 0; i < avgQualitiesSize; i++)
        sumQualities += avgQualities[i] * readLens[i];  // sum the qualities normalized by their read length

    return (double) sumQualities/getTotReadLen();
}

void InReads::report() {

    if (totReads > 0) {
        
        readLens.sort();
        
        std::cout << std::fixed; // disables scientific notation
        std::cout << std::setprecision(2); // 2 decimal points
        
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
        std::cout<<output("GC content %")<<gfa_round(computeGCcontent())<<"\n";
        std::cout<<output("Base composition (A:C:T:G)")<<totA<<":"<<totC<<":"<<totT<<":"<<totG<<"\n";
        std::cout<<output("Average read quality")<<getAvgQuality()<<"\n";
    }
}

void InReads::printReadLengths() {
    
    std::cout << std::fixed; // disables scientific notation
    std::cout << std::setprecision(2); // 2 decimal points

    if (userInput.sizeOutType == 's' || userInput.sizeOutType == 'h' || userInput.sizeOutType == 'c') {
        readLens.sort();
    }else{
        
        auto readLensTmp = readLens.all();
        
        for (auto read : readLensTmp)
            std::cout << read.first << "\n";
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
}

void InReads::printQualities() {
    
    std::cout << std::fixed; // disables scientific notation
    std::cout << std::setprecision(2); // 2 decimal points
    
    char qualityOut = userInput.qualityOut;
    
    auto readLensTmp = readLens.all();

    if (qualityOut == 'c'){
        for (auto read : readLensTmp)
            std::cout << read.second << "\n";
    }
    else if (qualityOut == 'l') { // l prints read lengths and qualities
        for (auto read : readLensTmp)
            std::cout << read.first << "," << read.second << "\n";
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
    
    if (streamOutput) {
        
        std::string ext = getFileExt(outputStream.file); // variable to handle output path and extension
        
        const static phmap::flat_hash_map<std::string,int> string_to_case{
            {"fasta",1},
            {"fa",1},
            {"fasta.gz",1},
            {"fa.gz",1},
            {"fastq",2},
            {"fq",2},
            {"fastq.gz",2},
            {"fq.gz",2}
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

void InReads::printTableCompressedBinary(std::string outFile) {
    
    std::ofstream ofs(outFile, std::fstream::trunc | std::ios::out | std::ios::binary);
    
    // write summary statistics
    ofs.write(reinterpret_cast<const char*>(&totA), sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&totC), sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&totG), sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&totT), sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&totN), sizeof(uint64_t));
    
    std::vector<std::pair<uint8_t,float>> &readLens8 = readLens.getReadLens8();
    std::vector<std::pair<uint16_t,float>> &readLens16 = readLens.getReadLens16();
    std::vector<std::pair<uint64_t,float>> &readLens64 = readLens.getReadLens64();
    
    // write vector lengths
    uint64_t len8 = readLens8.size(), len16 = readLens16.size(), len64 = readLens64.size();
    ofs.write(reinterpret_cast<const char*>(&len8), sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&len16), sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&len64), sizeof(uint64_t));
 
    // write vectors
    ofs.write(reinterpret_cast<const char*>(&readLens8[0]), len8 * sizeof(readLens8[0]));
    ofs.write(reinterpret_cast<const char*>(&readLens16[0]), len16 * sizeof(readLens16[0]));
    ofs.write(reinterpret_cast<const char*>(&readLens64[0]), len64 * sizeof(readLens64[0]));
    ofs.write(reinterpret_cast<const char*>(&avgQualities[0]), (len8 + len16 + len64) * sizeof(double));
    ofs.close();
}

void InReads::readTableCompressedBinary(std::string inFile) {
    
    // read
    std::ifstream ifs;
    ifs.open(inFile, std::ifstream::binary);
    
    // read summary statistics
    uint64_t A, C, G, T, N;
    ifs.read(reinterpret_cast<char*> (&A), sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*> (&C), sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*> (&G), sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*> (&T), sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*> (&N), sizeof(uint64_t));
    totA += A;
    totC += C;
    totG += G;
    totT += T;
    totN += N;
    
    uint64_t len8, len16, len64;
    ifs.read(reinterpret_cast<char*> (&len8), sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*> (&len16), sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*> (&len64), sizeof(uint64_t));
    
    // tmp vectors
    LenVector<float> readLensTmp;
    std::vector<float> avgQualitiesTmp;
    
    std::vector<std::pair<uint8_t,float>> &readLensTmp8 = readLensTmp.getReadLens8();
    std::vector<std::pair<uint16_t,float>> &readLensTmp16 = readLensTmp.getReadLens16();
    std::vector<std::pair<uint64_t,float>> &readLensTmp64 = readLensTmp.getReadLens64();
    
    readLensTmp8.resize(len8);
    readLensTmp16.resize(len16);
    readLensTmp64.resize(len64);
    avgQualitiesTmp.resize(len8 + len16 + len64);
    
    ifs.read(reinterpret_cast<char*> (&readLensTmp8[0]), len8 * sizeof(readLensTmp8[0]));
    ifs.read(reinterpret_cast<char*> (&readLensTmp16[0]), len16 * sizeof(readLensTmp16[0]));
    ifs.read(reinterpret_cast<char*> (&readLensTmp64[0]), len64 * sizeof(readLensTmp64[0]));
    ifs.read(reinterpret_cast<char*> (&avgQualitiesTmp[0]), (len8 + len16 + len64) * sizeof(double));
    
    // add to vector
    readLens.insert(readLensTmp);
    avgQualities.insert(std::end(avgQualities), std::begin(avgQualitiesTmp), std::end(avgQualitiesTmp));
    
    totReads += len8 + len16 + len64;
}
