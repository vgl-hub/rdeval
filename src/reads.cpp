#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>

#include <htslib/sam.h>
#include <htslib/thread_pool.h>

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

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )

void InRead::set(Log* threadLog, uint32_t uId, uint32_t iId, std::string seqHeader, std::string* seqComment, std::string* sequence, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, uint32_t seqPos, std::string* sequenceQuality, float avgQuality, std::vector<Tag>* inSequenceTags, uint64_t* N) {
    
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
    
    if (sequence != NULL && *sequence != "*") {
        this->setInSequence(sequence);
        threadLog->add("Segment sequence set");
    }
        
    if (sequenceQuality != NULL) {
        this->setInSequenceQuality(sequenceQuality);
        threadLog->add("Segment sequence quality set");
    }
    this->avgQuality = avgQuality;
    
    this->setACGT(A, C, G, T, N);
    threadLog->add("Increased ACGT counts");
    this->setLowerCount(lowerCount);
    threadLog->add("Increased total count of lower bases");
    
    if (sequence != NULL && *sequence != "*") {
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
    md5s.reserve(numFiles); // to avoid invalidating the vector during thread concurrency
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
        {"cram",2},
        {"rd",3}
    };
    
    for (uint32_t i = 0; i < numFiles; i++) {
        
        std::string file = userInput.file('r', i);
        std::string ext = getFileExt(file);
        if (ext != "rd") {
            md5s.push_back(std::make_pair(getFileName(file),std::string()));
            threadPool.queueJob([=]{ return computeMd5(file, md5s.back().second); });
        }
        
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
                                
                                size_t spacePos = newLine.find(" ");
                                seqHeader = newLine.substr(0, spacePos);
                                if (spacePos != std::string::npos)
                                    seqComment = newLine.substr(spacePos + 1);
                                else
                                    seqComment.clear();
                                
                                std::string* inSequence = new std::string;
								if (!getline(*stream, *inSequence, '>')) {
									fprintf(stderr, "Record appears truncated (%s). Exiting.\n", seqHeader.c_str());
									exit(EXIT_FAILURE);
								}
                                readBatch->sequences.push_back(new Sequence {seqHeader, seqComment, inSequence});
                                seqPos++;
                                processedLength += inSequence->size();
                                
                                if (processedLength > batchSize) {
                                    readBatch->batchN = batchN++;
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
                                size_t spacePos = newLine.find(" ");
                                seqHeader = newLine.substr(0, spacePos);
                                if (spacePos != std::string::npos)
                                    seqComment = newLine.substr(spacePos + 1);
                                else
                                    seqComment.clear();
                                
                                std::string* inSequence = new std::string;
								if (!getline(*stream, *inSequence)) {
									fprintf(stderr, "Record appears truncated (%s). Exiting.\n", seqHeader.c_str());
									exit(EXIT_FAILURE);
								}
								if (!getline(*stream, newLine)) {
									fprintf(stderr, "Record appears truncated (%s). Exiting.\n", seqHeader.c_str());
									exit(EXIT_FAILURE);
								}
                                std::string* inSequenceQuality = new std::string;
								if (!getline(*stream, *inSequenceQuality)) {
									fprintf(stderr, "Record appears truncated (%s). Exiting.\n", seqHeader.c_str());
									exit(EXIT_FAILURE);
								}
                                
                                readBatch->sequences.push_back(new Sequence {seqHeader, seqComment, inSequence, inSequenceQuality});
                                ++seqPos;
                                processedLength += inSequence->size();
                                
                                if (processedLength > batchSize) {
                                    readBatch->batchN = batchN++;
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
                    readBatch->batchN = batchN++; // process residual reads
                    lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
                    appendReads(readBatch);
                }
                break;
            }
            case 2: { // bam, cram
                
                Sequences* readBatch = new Sequences;
                samFile *fp_in = hts_open(userInput.file('r', i).c_str(),"r"); //open bam file
                bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
                bam1_t *bamdata = bam_init1(); //initialize an alignment
                
                tpool_read = {NULL, 0};
                tpool_read.pool = hts_tpool_init(userInput.decompression_threads);
                if (tpool_read.pool)
                    hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &tpool_read);
                else
                    lg.verbose("Failed to generate decompression threadpool with " + std::to_string(userInput.decompression_threads) + " threads. Continuing single-threaded");
                
                while(sam_read1(fp_in,bamHdr,bamdata) > 0) {
                    
                    uint32_t len = bamdata->core.l_qseq; // length of the read.
                    uint8_t *seq = bam_get_seq(bamdata); // seq string
                    std::string* inSequenceQuality = NULL;
                    
                    std::string* inSequence = new std::string;
                    inSequence->resize(len);
                    for(uint32_t i=0; i<len; ++i)
                        inSequence->at(i) = seq_nt16_str[bam_seqi(seq,i)]; //gets nucleotide id and converts them into IUPAC id.
                    
                    if (bam_get_qual(bamdata)[0] != (uint8_t)-1) {
                        inSequenceQuality = new std::string((char*)bam_get_qual(bamdata),len);
                        
                        for(uint32_t i=0; i<len; ++i)
                            inSequenceQuality->at(i) += 33;
                    }else{
                        uint8_t* tag = bam_aux_get(bamdata, "mq");
                        if (tag != 0)
                            inSequenceQuality = new std::string((char*)tag++,len);
                    }
                    readBatch->sequences.push_back(new Sequence {bam_get_qname(bamdata), std::string(), inSequence, inSequenceQuality});
                    seqPos++;
                    processedLength += inSequence->size();
                    
                    if (processedLength > batchSize) {
                        readBatch->batchN = batchN++;
                        lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
                        appendReads(readBatch);
                        readBatch = new Sequences;
                        processedLength = 0;

                    }
                    lg.verbose("Individual fastq sequence read: " + seqHeader);

                }
                readBatch->batchN = batchN++;
                lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
                appendReads(readBatch);
                bam_destroy1(bamdata);
                sam_close(fp_in);
                if (tpool_read.pool)
                    hts_tpool_destroy(tpool_read.pool);
                break;
            }
            case 3: { // rd
                    readTableCompressed(userInput.inFiles[i]);
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
}

float InReads::computeAvgQuality(std::string &sequenceQuality) {
    double sumQuality = 0;
    for (char &quality : sequenceQuality)
        sumQuality += pow(10, -(double(quality) - 33)/10);
    return (float) -10 * std::log10(sumQuality/sequenceQuality.size());
}

void InReads::initDictionaries() {
    
    std::ifstream includeFile(userInput.inBedInclude);
    std::string key;

    while (includeFile >> key)
        includeList.insert(key);

    includeFile.close();
    
    std::ifstream excludeFile(userInput.inBedExclude);

    while (excludeFile >> key)
        excludeList.insert(key);

    excludeFile.close();
}

void InReads::initFilters() {
    
    bool cannotParse = false;
    
    std::istringstream stream(userInput.filter); // get l,q and logical operator in between
    std::string token1;
    std::string token2;
    std::string token3;
    stream >> token1 >> token2 >> token3;

    if ((userInput.filter.find('l') == std::string::npos && userInput.filter.find('q') == std::string::npos) || // no filter found
        (token1.size() < 3) || // first filter too short
        (token2.size() != 0 && token3.size() == 0) || // second filter missing
        (token2.size() == 0 && token3.size() != 0) || // missing operator
        (token2.size() > 1) // malformatted operator
    )
        cannotParse = true;
    
    if (cannotParse){
        fprintf(stderr, "Could not parse filter. Terminating.\n");
        exit(EXIT_FAILURE);
    }
    
    std::string lFilter, qFilter;
    if (token2.size())
        logicalOperator = token2[0];
    
    if (token1[0] == 'l') {
        lFilter = token1;
        lSign = lFilter[1];
        l = stoi(lFilter.substr(2));
    } else if (token1[0] == 'q') {
        qFilter = token1;
        qSign = qFilter[1];
        q = stoi(qFilter.substr(2));
    }
    
    if (token3.size() != 0) {
        if (token3[0] == 'l') {
            lFilter = token3;
            lSign = lFilter[1];
            l = stoi(lFilter.substr(2));
        }else if (token3[0] == 'q') {
            qFilter = token3;
            qSign = qFilter[1];
            q = stoi(qFilter.substr(2));
        }
    }
}

void InReads::filterRecords() {
	
	uint64_t readCount = readLens.size();
	LenVector<float> readLensFiltered;
	for (uint64_t i = 0; i < readCount; ++i) {
		if (!applyFilter(readLens[i].first, readLens[i].second))
			readLensFiltered.push_back(std::make_pair(readLens[i].first, readLens[i].second));
	}
	if (readLens.size() != readLensFiltered.size()) { // if the read counts after filtering is different these counts are invalidated
		totA = 0;
		totC = 0;
		totG = 0;
		totT = 0;
		totN = 0;
	}
	readLens = readLensFiltered; // overwrite the counts
	totReads = readLens.size();
}



inline bool InReads::filterRead(Sequence* sequence) {
    
    uint64_t size = sequence->sequence->size();
    float avgQuality = 0;
    if (sequence->sequenceQuality != NULL)
        avgQuality = computeAvgQuality(*sequence->sequenceQuality);
    
    return applyFilter(size, avgQuality);
}

inline bool InReads::applyFilter(uint64_t size, float avgQuality) {
    
    bool lFilter = ((lSign == '0') ||
                    ((lSign == '>') && (size > l)) ||
                    ((lSign == '<') && (size < l)) ||
                    ((lSign == '=') && (size == l))
                    );
    
    bool qFilter = ((qSign == '0') ||
                    ((qSign == '>') && (avgQuality > q)) ||
                    ((qSign == '<') && (avgQuality < q)) ||
                    ((qSign == '=') && (avgQuality == q))
                    );
    
    if ((logicalOperator == '0' && (lSign != '0' && lFilter)) ||
        (logicalOperator == '0' && (qSign != '0' && qFilter)) ||
        (logicalOperator == '|' && (lFilter || qFilter)) ||
        (logicalOperator == '&' && (lFilter && qFilter))
       )
        return false;
    return true;
}

float newRand() {
    return (static_cast <float> (random()) / static_cast <float> (RAND_MAX));
}

bool InReads::traverseInReads(Sequences* readBatch) { // traverse the read

    Log threadLog;
    threadLog.setId(readBatch->batchN);
    std::vector<bam1_t*> inReadsBatch;
    std::vector<InRead*> inReadsSummaryBatch;
    uint32_t readN = 0;
    LenVector<float> readLensBatch;
    InRead* read;

    uint64_t batchA = 0, batchT = 0, batchC = 0, batchG = 0, batchN = 0;
    bool include = includeList.size();
    bool exclude = excludeList.size();
    bool filter = userInput.filter != "none" ? true : false;
    bool sample = userInput.ratio < 1 ? true : false; // read subsampling

    for (Sequence* sequence : readBatch->sequences) {
        
        if (include) {
            if (includeList.find(sequence->header) == includeList.end())
                continue;
        }
        if (exclude) {
            if (excludeList.find(sequence->header) != excludeList.end())
                continue;
        }
        if (filter) {
            bool filtered = filterRead(sequence);
            if (filtered)
                continue;
        }
        if (sample) {
            if (newRand() > userInput.ratio)
                continue;
        }
        read = traverseInRead(&threadLog, sequence, readBatch->batchN+readN++);
        std::pair<uint64_t, float> lenQual(read->getA()+read->getC()+read->getG()+read->getT()+read->getN(), read->avgQuality);
        readLensBatch.push_back(lenQual);
        batchA += read->getA();
        batchT += read->getT();
        batchC += read->getC();
        batchG += read->getG();
        batchN += read->getN();

        if (streamOutput) {
            
            bam1_t *q;
            if (!(q = bam_init1())) {
                printf("Failed to initialize bamdata\n");
                exit(EXIT_FAILURE);
            }

            if (read->inSequenceQuality == NULL)
                read->inSequenceQuality = new std::string(read->inSequence->size(),'!');
            
            if (bam_set1(q, strlen(read->seqHeader.c_str()), read->seqHeader.c_str(), BAM_FUNMAP, -1, -1, 0, 0, NULL, -1, -1, 0, strlen(read->inSequence->c_str()), read->inSequence->c_str(), NULL, 0) < 0) {
                printf("Failed to set data\n");
                exit(EXIT_FAILURE);
            }
            uint8_t *s = bam_get_seq(q);
            s = bam_get_qual(q);
            for (size_t i = 0; i < read->inSequenceQuality->size(); ++i)
                s[i] = read->inSequenceQuality->at(i) - 33;
            
            inReadsBatch.push_back(q);
        }
		if (userInput.content_flag)
			inReadsSummaryBatch.push_back(read);
		else
			delete read;
    }
    
    {
        std::unique_lock<std::mutex> lck(mtx);
        readBatches.emplace_back(inReadsBatch,readBatch->batchN);
        readSummaryBatches.emplace_back(inReadsSummaryBatch,readBatch->batchN);
        delete readBatch;
        readLens.insert(readLensBatch);
        totA+=batchA;
        totT+=batchT;
        totC+=batchC;
        totG+=batchG;
        totN+=batchN;
        totReads += readN;
        logs.push_back(threadLog);
    }
    writerMutexCondition.notify_one();
    return true;
}

InRead* InReads::traverseInRead(Log* threadLog, Sequence* sequence, uint32_t seqPos) { // traverse a single read

    std::vector<std::pair<uint64_t, uint64_t>> bedCoords;
    if(userInput.hc_cutoff != -1) {
        homopolymerCompress(sequence->sequence, bedCoords, userInput.hc_cutoff);
        delete sequence->sequenceQuality; // sequence quality not meaningful when compressed
        sequence->sequenceQuality = NULL;
    }
    uint64_t A = 0, C = 0, G = 0, T = 0, N = 0, lowerCount = 0;
    float avgQuality = 0;
    
    for (char &base : *sequence->sequence) {
        
        if (islower(base))
            ++lowerCount;
                
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
    if (sequence->sequenceQuality != NULL)
        avgQuality = computeAvgQuality(*sequence->sequenceQuality);

    // operations on the segment
    InRead* inRead = new InRead;
    
    if (userInput.content_flag && !streamOutput) {
        delete sequence->sequence;
        sequence->sequence = NULL;
        if (sequence->sequenceQuality != NULL) {
            delete sequence->sequence;
            sequence->sequenceQuality = NULL;
        }
    }
    inRead->set(threadLog, 0, 0, sequence->header, &sequence->comment, sequence->sequence, &A, &C, &G, &T, &lowerCount, seqPos, sequence->sequenceQuality, avgQuality, NULL, &N);
    sequence->sequence = NULL;
    sequence->sequenceQuality = NULL;
    return inRead;
}

uint64_t InReads::getTotReadLen() {
	
	uint64_t sum = totA + totC + totG + totT + totN;
	if (!sum) {
		uint64_t readCount = readLens.size();
		for (uint64_t i = 0; i < readCount; ++i)
			sum += readLens[i].first;
	}
    return sum;
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

void InReads::evalNstars() { // clean up once len-vector iterator is available, still expensive

    uint64_t sum = 0, totLen = getTotReadLen();
    uint8_t N = 1;
    
    for(unsigned int i = 0; i < readLens.size(); ++i) { // for each length
        sum += readLens[i].first; // increase sum
        while (sum >= ((double) totLen / 10 * N) && N<= 10) { // conditionally add length.at or pos to each N/L* bin
            
            readNstars[N-1] = readLens[i].first;
            readLstars[N-1] = i + 1;
            N += 1;
        }
    }
}

uint64_t InReads::getSmallestRead() {
    return readLens.back();
}

uint64_t InReads::getLargestRead() {
    return readLens.front();
}

double InReads::getAvgQuality(){

    double sumQualities = 0, avgQualitiesSize=readLens.size();

    for (uint64_t i = 0; i < avgQualitiesSize; ++i)
        sumQualities += readLens[i].first * pow(10, -readLens[i].second/10);  // sum the qualities normalized by their read length

    return -10 * std::log10(sumQualities/getTotReadLen());
}

void InReads::report() {

    if (totReads > 0) {
		
		if (userInput.filter != "none" && getFileExt(userInput.file('r', 0)) == "rd")
			filterRecords();
        
        readLens.sort();
        
        std::cout << std::fixed; // disables scientific notation
        std::cout << std::setprecision(2); // 2 decimal points
        
        if (!tabular_flag)
            std::cout<<output("+++Read summary+++")<<"\n";
        std::cout<<output("# reads")<<totReads<<std::endl;
        std::cout<<output("Total read length")<<getTotReadLen()<<std::endl;
        std::cout<<output("Average read length") << gfa_round(computeAvgReadLen())<<std::endl;
        evalNstars(); // read N* statistics
        std::cout<<output("Read N50")<<getReadN50()<<std::endl;
        std::cout<<output("Smallest read length")<<getSmallestRead()<<std::endl;
        std::cout<<output("Largest read length")<<getLargestRead()<<std::endl;
        std::cout<<output("Coverage")<<gfa_round((double)getTotReadLen()/userInput.gSize)<<std::endl;
        std::cout<<output("GC content %")<<gfa_round(computeGCcontent())<<std::endl;
        std::cout<<output("Base composition (A:C:T:G)")<<totA<<":"<<totC<<":"<<totT<<":"<<totG<<std::endl;
        std::cout<<output("Average per base quality")<<std::abs(getAvgQuality())<<std::endl;
    }
}

void InReads::printReadLengths() {
    
    std::cout << std::fixed; // disables scientific notation
    std::cout << std::setprecision(2); // 2 decimal points
    
    uint64_t readCount = readLens.size();
    bool noFilter = userInput.filter == "none" ? true : false;
    if(userInput.sizeOutType == 's')
        readLens.sort();
    
    if (userInput.sizeOutType == 'u' || userInput.sizeOutType == 's') {
        
        for (uint64_t i = 0; i < readCount; ++i) {
            if (noFilter || !applyFilter(readLens[i].first, readLens[i].second))
                std::cout << readLens[i].first << "\n";
        }
        
    }else if(userInput.sizeOutType == 'h') {
        
        phmap::parallel_flat_hash_map<uint64_t, uint64_t> hist;
        
        for (uint64_t i = 0; i < readCount; ++i) {
            if (noFilter || !applyFilter(readLens[i].first, readLens[i].second))
                ++hist[readLens[i].first];
        }
        std::vector<std::pair<uint64_t, uint64_t>> table(hist.begin(), hist.end()); // converts the hashmap to a table
        std::sort(table.begin(), table.end());
        
        for (auto pair : table)
            std::cout<<pair.first<<"\t"<<pair.second<<"\n";
        
    }else if(userInput.sizeOutType == 'c') {

        phmap::parallel_flat_hash_map<uint64_t, uint64_t> hist;
        
        for (uint64_t i = 0; i < readCount; ++i) {
            if (noFilter || !applyFilter(readLens[i].first, readLens[i].second))
                ++hist[readLens[i].first];
        }
        std::vector<std::pair<uint64_t, uint64_t>> table(hist.begin(), hist.end());
        std::sort(table.begin(), table.end());
        uint64_t totReadLen = getTotReadLen(), sum = 0;
        
        for (auto pair : table) {
            std::cout<<+pair.first<<"\t"<<+pair.second<<"\t"<<+pair.first*pair.second<<"\t"<<+(totReadLen - sum)<<"\n";
            sum += pair.first*pair.second;
        }
    }
}

void InReads::printQualities() {
    
    std::cout << std::fixed; // disables scientific notation
    std::cout << std::setprecision(2); // 2 decimal points
    
    uint64_t readCount = readLens.size();
    
    bool noFilter = userInput.filter == "none" ? true : false;
    char qualityOut = userInput.qualityOut;

    if (qualityOut == 'q'){
        for (uint64_t i = 0; i < readCount; ++i) {
            if (noFilter || !applyFilter(readLens[i].first, readLens[i].second))
                std::cout << readLens[i].second << "\n";
        }
    }
    else if (qualityOut == 'a') { // a prints read lengths and qualities
        for (uint64_t i = 0; i < readCount; ++i) {
            if (noFilter || !applyFilter(readLens[i].first, readLens[i].second))
                std::cout << readLens[i].first << "," << readLens[i].second << "\n";
        }
    }
}
void InReads::printContent() {
    
    std::cout<<"Header\tComment\tLength\tA\tC\tG\tT\tN\tGC\tAverage Quality\n";
    for (std::pair<std::vector<InRead*>,uint32_t> inReads : readSummaryBatches){
        for (InRead* read : inReads.first){
            uint64_t A = read->getA(), C = read->getC(), G = read->getG(), T = read->getT(), N = read->getN(), total = A+C+G+T+N;
            std::cout<<read->seqHeader<<"\t"<<read->seqComment<<"\t"<<total<<"\t"<<A<<"\t"<<C<<"\t"<<G<<"\t"<<T<<"\t"<<N<<"\t"<<gfa_round((float)(G+C)/total)<<"\t"<<read->avgQuality<<"\n";
            delete read;
        }
    }
}

void dump_read(bam1_t* b) {
    printf("->core.tid:(%d)\n", b->core.tid);
    printf("->core.pos:(%ld)\n", static_cast<unsigned long>(b->core.pos));
    printf("->core.bin:(%d)\n", b->core.bin);
    printf("->core.qual:(%d)\n", b->core.qual);
    printf("->core.l_qname:(%d)\n", b->core.l_qname);
    printf("->core.flag:(%d)\n", b->core.flag);
    printf("->core.n_cigar:(%d)\n", b->core.n_cigar);
    printf("->core.l_qseq:(%d)\n", b->core.l_qseq);
    printf("->core.mtid:(%d)\n", b->core.mtid);
    printf("->core.mpos:(%ld)\n", static_cast<unsigned long>(b->core.mpos));
    printf("->core.isize:(%ld)\n", static_cast<unsigned long>(b->core.isize));
    if (b->data) {
        printf("->data:");
        int i;
        for (i = 0; i < b->l_data; ++i) {
            printf("%x ", b->data[i]);
        }
        printf("\n");
    }
    if (b->core.l_qname) {
        printf("qname: %s\n",bam_get_qname(b));
    }
    if (b->core.l_qseq) {
        printf("qseq:");
        int i;
        for (i = 0; i < b->core.l_qseq; ++i) {
            printf("%c", seq_nt16_str[bam_seqi(bam_get_seq(b),i)]);
        }
        printf("\n");
        printf("qual:");
        uint8_t *s = bam_get_qual(b);
        for (i = 0; i < b->core.l_qseq; ++i) {
            printf("%c",s[i] + 33);
        }
        printf("\n");

    }

    if (bam_get_l_aux(b)) {
        uint32_t i = 0;
        uint8_t* aux = bam_get_aux(b);

        while (i < bam_get_l_aux(b)) {
            printf("%.2s:%c:",aux+i,*(aux+i+2));
            i += 2;
            switch (*(aux+i)) {
                case 'Z':
                    while (*(aux+1+i) != '\0') { putc(*(aux+1+i), stdout); ++i; }
                    break;
            }
            putc('\n',stdout);
            ++i;++i;
        }
    }
    printf("\n");
}

void InReads::initStream() {
    
    if (streamOutput) {
        
        if(userInput.outFiles.size() > 1) {
            fprintf(stderr, "Error: rdeval does not support more than one output at a time. Terminating.\n");
            exit(EXIT_FAILURE);
        }
        
        std::string ext = getFileExt(userInput.outFiles[0]); // variable to handle output path and extension
        const static phmap::flat_hash_map<std::string,int> string_to_case{
            {"fasta",1},
            {"fa",1},
            {"fasta.gz",1},
            {"fa.gz",1},
            {"fastq",2},
            {"fq",2},
            {"fastq.gz",2},
            {"fq.gz",2},
            {"bam",3},
            {"cram",4}
        };
        
        switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
            case 1:   // fasta[.gz]
            case 2: { // fastq[.gz]
                char mode[4] = "w";
                if (sam_open_mode(mode + 1, userInput.outFiles[0].c_str(), NULL) < 0) {
                    printf("Invalid file name\n");
                    exit(EXIT_FAILURE);
                }
                if (!(fp = sam_open(userInput.outFiles[0].c_str(), mode))) {
                    printf("Could not open %s\n", userInput.outFiles[0].c_str());
                    exit(EXIT_FAILURE);
                }
                break;
            }
            case 3: {  // bam
                fp = sam_open(userInput.outFiles[0].c_str(),"wb");
                break;
            }
            case 4: {  // cram
                htsFormat fmt4 = {sequence_data, cram, {3, 1}, gzip, 6, NULL};
                hts_parse_format(&fmt4, "cram,no_ref=1");
                fp = sam_open_format(userInput.outFiles[0].c_str(), "wc", &fmt4);
                break;
            }
        }
        if (bam)
            writeHeader();
        writer = std::thread(&InReads::writeToStream, this);
    }
}

void InReads::closeStream() {
    
    if (writer.joinable()) {
        {
            std::unique_lock<std::mutex> lck(mtx);
            streamOutput = false;
        }
        writerMutexCondition.notify_one();
        writer.join();
    }
}

void InReads::writeHeader() {
    
    const char init_header[] = "@HD\tVN:1.4\tSO:unknown\n";
    hdr = bam_hdr_init();
    hdr->l_text = strlen(init_header);
    hdr->text = strdup(init_header);
    hdr->n_targets = 0;
    std::ignore = sam_hdr_write(fp,hdr);
}

void InReads::closeBam() {
    bam_hdr_destroy(hdr);
}

void InReads::writeToStream() {

    uint64_t batchCounter = 0;
    std::vector<std::pair<std::vector<bam1_t*>,uint32_t>> readBatchesCpy;
    
    tpool_write = {NULL, 0}; // init htslib threadpool
    tpool_write.pool = hts_tpool_init(userInput.compression_threads);
    if (tpool_write.pool)
        hts_set_opt(fp, HTS_OPT_THREAD_POOL, &tpool_write);
    else
        lg.verbose("Failed to generate compression threadpool with " + std::to_string(userInput.compression_threads) + " threads. Continuing single-threaded");
    
    while (true) {
        
        {
            std::unique_lock<std::mutex> lck(mtx);
            writerMutexCondition.wait(lck, [this, &readBatchesCpy] {
                return !streamOutput || readBatches.size() || readBatchesCpy.size();
            });
            if (!streamOutput && !readBatches.size() && !readBatchesCpy.size()) {
                if (bam)
                    closeBam();
                sam_close(fp); // close file
                if (tpool_write.pool)
                    hts_tpool_destroy(tpool_write.pool);
                return;
            }
            readBatchesCpy.insert(readBatchesCpy.begin(), readBatches.begin(), readBatches.end());
            readBatches.clear();
        }
                
        for (std::vector<std::pair<std::vector<bam1_t*>,uint32_t>>::iterator it = readBatchesCpy.begin(); it != readBatchesCpy.end();) {
            
            auto &inReads = *it;
            
            if (inReads.second != batchCounter) {
                ++it;
                continue;
            }
            lg.verbose("Writing read batch " + std::to_string(inReads.second) + " to file (" + std::to_string(inReads.first.size())  + ")");
                
            for (bam1_t* q : inReads.first){
                if (sam_write1(fp, hdr, q) < 0) {
                    printf("Failed to write data\n");
                    exit(EXIT_FAILURE);
                }
                bam_destroy1(q);
            }
            it = readBatchesCpy.erase(it);
            ++batchCounter;
        }
    }
}

void InReads::printTableCompressed(std::string outFile) {
	
    // compute buffer size
    std::vector<std::pair<uint8_t,float>> &readLens8 = readLens.getReadLens8();
    std::vector<std::pair<uint16_t,float>> &readLens16 = readLens.getReadLens16();
    std::vector<std::pair<uint64_t,float>> &readLens64 = readLens.getReadLens64();
    uint64_t len8 = readLens8.size(), len16 = readLens16.size(), len64 = readLens64.size();
    
    uLong sourceLen = sizeof(uint64_t) * (5 + 3) + len8 * sizeof(readLens8[0]) + len16 * sizeof(readLens16[0]) + len64 * sizeof(readLens64[0]); // ACGTN + len8,len16,len64 + vectors
    Bytef *source = new Bytef[sourceLen];
    uLong destLen = compressBound(sourceLen);
    Bytef *dest = new Bytef[destLen];

    unsigned char* ptr = source;

    // write ACGTN counts
    memcpy(ptr, &totA, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &totC, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &totG, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &totT, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &totN, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    
    // write vector lengths
    memcpy(ptr, &len8, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &len16, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &len64, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
 
    // write vectors
    memcpy(ptr, &readLens8[0], len8 * sizeof(readLens8[0]));
    ptr += len8 * sizeof(readLens8[0]);
    memcpy(ptr, &readLens16[0], len16 * sizeof(readLens16[0]));
    ptr += len16 * sizeof(readLens16[0]);
    memcpy(ptr, &readLens64[0], len64 * sizeof(readLens64[0]));
    ptr += len64 * sizeof(readLens64[0]);

    compress(dest, &destLen, source, sourceLen);
    delete[] source;
    
    std::ofstream ofs(outFile, std::fstream::trunc | std::ios::out | std::ios::binary);
    
    // write md5s
    uint32_t md5sN = md5s.size();
    uint16_t stringSize;
    ofs.write(reinterpret_cast<const char*>(&md5sN), sizeof(uint32_t)); // <-- these could presumably be also gzipped (need to work on the R interface too)
    for (auto md5 : md5s) {
        stringSize = md5.first.size();
        ofs.write(reinterpret_cast<const char*>(&stringSize), sizeof(uint16_t));
        ofs.write(reinterpret_cast<const char*>(md5.first.c_str()), sizeof(char) * stringSize);
        stringSize = md5.second.size();
        ofs.write(reinterpret_cast<const char*>(&stringSize), sizeof(uint16_t));
        ofs.write(reinterpret_cast<const char*>(md5.second.c_str()), sizeof(char) * stringSize);
    }
    ofs.write(reinterpret_cast<const char*>(&sourceLen), sizeof(uint64_t)); // output the decompressed file size
    ofs.write(reinterpret_cast<const char*>(dest), destLen * sizeof(Bytef)); // output compressed data
    ofs.close();
    delete[] dest;
}

void InReads::readTableCompressed(std::string inFile) {
	
	if (userInput.filter != "none" && getFileExt(userInput.file('r', 0)) == "rd")
		filterRecords();
    
    // read
    std::ifstream ifs(inFile, std::ios::binary | std::ios::ate); // compute file size
    std::streamsize fileSize = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    
    uint32_t md5sN;
    uint16_t stringSize;
    ifs.read(reinterpret_cast<char*>(&md5sN), sizeof(uint32_t));
    
    for (uint32_t i = 0; i < md5sN; ++i) {
        std::string filename, md5;
        ifs.read(reinterpret_cast<char*>(&stringSize), sizeof(uint16_t));
        filename.resize(stringSize);
        ifs.read(reinterpret_cast<char*>(&filename[0]), sizeof(char) * stringSize);
        ifs.read(reinterpret_cast<char*>(&stringSize), sizeof(uint16_t));
        md5.resize(stringSize);
        ifs.read(reinterpret_cast<char*>(&md5[0]), sizeof(char) * stringSize);
        md5s.push_back(std::make_pair(filename,md5));
    }
    
    uLongf decompressedSize;
    ifs.read(reinterpret_cast<char*> (&decompressedSize), sizeof(uint64_t)); // read gz-uncompressed size
    Bytef *data = new Bytef[decompressedSize];
    
    uLong compressedSize = fileSize - sizeof(uint64_t); // it seems this should subtract the header too, probably it should be fixed subtracting the md5s
    Bytef *gzData = new Bytef[compressedSize];
    ifs.read(reinterpret_cast<char*> (gzData), compressedSize * sizeof(Bytef)); // read gzipped data
    uncompress(data, &decompressedSize, gzData, compressedSize);
    delete[] gzData;
    
    // read summary statistics
    unsigned char* ptr = data;
    
    uint64_t A, C, G, T, N;
    memcpy(&A, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&C, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&G, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&T, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&N, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    totA += A;
    totC += C;
    totG += G;
    totT += T;
    totN += N;

    uint64_t len8, len16, len64;
    memcpy(&len8, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&len16, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&len64, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);

    // tmp vectors
    LenVector<float> readLensTmp;

    std::vector<std::pair<uint8_t,float>> &readLensTmp8 = readLensTmp.getReadLens8();
    std::vector<std::pair<uint16_t,float>> &readLensTmp16 = readLensTmp.getReadLens16();
    std::vector<std::pair<uint64_t,float>> &readLensTmp64 = readLensTmp.getReadLens64();

    readLensTmp8.resize(len8);
    readLensTmp16.resize(len16);
    readLensTmp64.resize(len64);
    // prefer this to memcpy as std::copy is safer
    std::pair<uint8_t,float>* p8 = reinterpret_cast<std::pair<uint8_t,float>*>(ptr);
    std::copy(p8, p8+len8, readLensTmp8.begin());
    ptr += len8 * sizeof(readLensTmp8[0]);
    std::pair<uint16_t,float>* p16 = reinterpret_cast<std::pair<uint16_t,float>*>(ptr);
    std::copy(p16, p16+len16, readLensTmp16.begin());
    ptr += len16 * sizeof(readLensTmp16[0]);
    std::pair<uint64_t,float>* p64 = reinterpret_cast<std::pair<uint64_t,float>*>(ptr);
    std::copy(p64, p64+len64, readLensTmp64.begin());
    ptr += len64 * sizeof(readLensTmp64[0]);
    
    delete[] data;
    
    // add to vector
    readLens.insert(readLensTmp);

    totReads += len8 + len16 + len64;
}

void InReads::printMd5() {
    for (auto md5 : md5s)
        std::cout<<md5.first<<": "<<md5.second<<std::endl;
}
