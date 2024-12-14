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

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )


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
        
        std::string file = userInput.file('r', i);
        std::string ext = getFileExt(file);
        if (ext != "rd") {
            std::string *md5 = new std::string;
            md5s.push_back(std::make_pair(getFileName(file),md5));
            threadPool.queueJob([=]{ return computeMd5(file, md5); });
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
                
                Sequences* readBatch = new Sequences;
                samFile *fp_in = hts_open(userInput.file('r', i).c_str(),"r"); //open bam file
                bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
                bam1_t *bamdata = bam_init1(); //initialize an alignment

//                char *tar = bamHdr->text;
//                printf("%s\n",tar);
                
                while(sam_read1(fp_in,bamHdr,bamdata) > 0) {
                    
                    uint32_t len = bamdata->core.l_qseq; //length of the read.
                    uint8_t *q = bam_get_seq(bamdata); //quality string
                    char *qseq = (char *)malloc(len);
                    char *qual = (char *)malloc(len);
                    
                    for(uint32_t i=0; i<len; ++i)
                        qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
                    
                    std::string* inSequence = new std::string(qseq, len);
                    free(qseq);
                    
                    for(int i=0; i<bamdata->core.l_qseq; ++i)
                        qual[i] = (char) bam_get_qual(bamdata)[i] + 33;
                    
                    std::string* inSequenceQuality = new std::string(qual, len);
                    free(qual);
                    
                    readBatch->sequences.push_back(new Sequence {bam_get_qname(bamdata), std::string(), inSequence, inSequenceQuality});
                    seqPos++;
                    processedLength += inSequence->size();
                    
                    if (processedLength > batchSize) {
                        readBatch->batchN = seqPos/batchSize;
                        lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
                        appendReads(readBatch);
                        readBatch = new Sequences;
                        processedLength = 0;

                    }
                    lg.verbose("Individual fastq sequence read: " + seqHeader);

                }
                readBatch->batchN = seqPos/batchSize + 1;
                lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));
                appendReads(readBatch);
                bam_destroy1(bamdata);
                sam_close(fp_in);
                break;
            }
            case 3: { // cram
                fprintf(stderr, "cram currently not supported.\n");
                exit(EXIT_FAILURE);
                break;
            }
            case 4: { // rd
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

        if (streamOutput)
            inReadsBatch.push_back(read);
        else
            delete read;
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    readBatches.emplace_back(inReadsBatch,readBatch->batchN);
    delete readBatch;
    readLens.insert(readLensBatch);
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

    uint64_t sumQualities = 0, avgQualitiesSize=readLens.size();

    for (uint64_t i = 0; i < avgQualitiesSize; i++)
        sumQualities += readLens[i].first * readLens[i].second;  // sum the qualities normalized by their read length

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
            if (readLens[i].first == readLens[i+1].first) {
                count += 1;
            }
            else if (readLens[i].first != readLens[i+1].first) {
                std::cout << readLens[i].first << "," << count << "\n";
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
            if (readLens[i].first == readLens[i+1].first) {
                count += 1;
            }
            else if (readLens[i].first != readLens[i+1].first) {
                sizexCount = (readLens[i].first * count);
                counts.push_back(count);
                sizexCounts.push_back(sizexCount);
                uniqReadLens.push_back(readLens[i].first);
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

void dump_read(bam1_t* b) {
    printf("->core.tid:(%d)\n", b->core.tid);
    printf("->core.pos:(%lld)\n", b->core.pos);
    printf("->core.bin:(%d)\n", b->core.bin);
    printf("->core.qual:(%d)\n", b->core.qual);
    printf("->core.l_qname:(%d)\n", b->core.l_qname);
    printf("->core.flag:(%d)\n", b->core.flag);
    printf("->core.n_cigar:(%d)\n", b->core.n_cigar);
    printf("->core.l_qseq:(%d)\n", b->core.l_qseq);
    printf("->core.mtid:(%d)\n", b->core.mtid);
    printf("->core.mpos:(%lld)\n", b->core.mpos);
    printf("->core.isize:(%lld)\n", b->core.isize);
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
            {"fq.gz",2},
            {"bam",3}
        };
        
        std::vector<std::pair<std::vector<InRead*>,uint32_t>> readBatchesCpy;
        {
            std::unique_lock<std::mutex> lck(mtx);
            readBatchesCpy = {readBatches.begin() + batchCounter-1, readBatches.end()};
        }
        
        switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
                
            case 1:  { // fasta[.gz]
                
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
            case 2:  { // fastq[.gz]
                
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
            case 3: { // bam
                
                const char init_header[] = "@HD\tVN:1.4\tSO:unknown\n";
                htsFile *fp = sam_open(outputStream.file.c_str(),"wb");

                // write header
                bam_hdr_t *hdr = bam_hdr_init();
                hdr->l_text = strlen(init_header);
                hdr->text = strdup(init_header);
                hdr->n_targets = 0;
                sam_hdr_write(fp,hdr);
                
                for (std::pair<std::vector<InRead*>,uint32_t> inReads : readBatchesCpy) {
                    
                    if (inReads.second > batchCounter)
                        continue;
                    
                    lg.verbose("Writing read batch " + std::to_string(inReads.second) + " to file (" + std::to_string(inReads.first.size())  + ")");
                        
                    for (InRead* read : inReads.first) { // main loop, iter through each fastq records
                        
                        std::string *quality;
                        if (read->inSequenceQuality != NULL)
                            quality = read->inSequenceQuality;
                        else
                            quality = new std::string('!', read->inSequence->size());
                        
                        bam1_t *q = bam_init1();
                        //`q->data` structure: qname-cigar-seq-qual-aux
                        q->l_data = read->seqHeader.size()+1+(int)(1.5*read->inSequence->size()+(read->inSequence->size() % 2 != 0)); // +1 includes the tailing '\0'
                        if (q->m_data < (uint32_t)q->l_data) {
                            q->m_data = q->l_data;
                            kroundup32(q->m_data);
                            q->data = (uint8_t*)realloc(q->data, q->m_data);
                        }
                        q->core.flag = BAM_FMUNMAP;
                        q->core.l_qname = read->seqHeader.size()+1; // +1 includes the tailing '\0'
                        q->core.l_qseq = read->inSequence->size();
                        q->core.n_cigar = 0; // we have no cigar sequence
                        // no flags for unaligned reads
                        q->core.tid = -1;
                        q->core.pos = -1;
                        q->core.mtid = -1;
                        q->core.mpos = -1;
                        memcpy(q->data, read->seqHeader.c_str(), q->core.l_qname); // first set qname
                        uint8_t *s = bam_get_seq(q);
                        for (int i = 0; i < q->core.l_qseq; ++i)
                            bam1_seq_seti(s, i, seq_nt16_table[read->inSequence->at(i)]);
                        s = bam_get_qual(q);
                        for (size_t i = 0; i < quality->size(); ++i)
                            s[i] = quality->at(i) - 33;
                        dump_read(q);
                        int ret = sam_write1(fp, hdr, q);
                        bam_destroy1(q);
                        if (read->inSequenceQuality != NULL)
                            delete quality;
                    }
                    ++batchCounter;
                }
                sam_close(fp); // close bam file
                break;
            }
        }
    }
}

void InReads::printTableCompressed(std::string outFile) {
    
    // compute buffer size
    std::vector<std::pair<uint8_t,float>> &readLens8 = readLens.getReadLens8();
    std::vector<std::pair<uint16_t,float>> &readLens16 = readLens.getReadLens16();
    std::vector<std::pair<uint64_t,float>> &readLens64 = readLens.getReadLens64();
    uint64_t len8 = readLens8.size(), len16 = readLens16.size(), len64 = readLens64.size();
    
    uLong sourceLen = sizeof(uint64_t) * (5 + 3) + len8 * sizeof(readLens8[0]) + len16 * sizeof(readLens16[0]) + sizeof(readLens64[0]) + (len8 + len16 + len64) * sizeof(float); // ACGTN + len8,len16,len64 + vectors + qualities
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
    ofs.write(reinterpret_cast<const char*>(&md5sN), sizeof(uint32_t));
    for (auto md5 : md5s) {
        stringSize = md5.first.size();
        ofs.write(reinterpret_cast<const char*>(&stringSize), sizeof(uint16_t));
        ofs.write(reinterpret_cast<const char*>(md5.first.c_str()), sizeof(char) * stringSize);
        stringSize = md5.second->size();
        ofs.write(reinterpret_cast<const char*>(&stringSize), sizeof(uint16_t));
        ofs.write(reinterpret_cast<const char*>(md5.second->c_str()), sizeof(char) * stringSize);
    }
    ofs.write(reinterpret_cast<const char*>(&sourceLen), sizeof(uint64_t)); // output the decompressed file size
    ofs.write(reinterpret_cast<const char*>(dest), destLen * sizeof(Bytef)); // output compressed data
    ofs.close();
    delete[] dest;
}

void InReads::readTableCompressed(std::string inFile) {
    
    // read
    std::ifstream ifs(inFile, std::ios::binary | std::ios::ate); // compute file size
    std::streamsize fileSize = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    
    uint32_t md5sN;
    uint16_t stringSize;
    ifs.read(reinterpret_cast<char*>(&md5sN), sizeof(uint32_t));
    
    for (uint32_t i = 0; i < md5sN; ++i) {
        std::string filename;
        std::string *md5 = new std::string;
        ifs.read(reinterpret_cast<char*>(&stringSize), sizeof(uint16_t));
        filename.resize(stringSize);
        ifs.read(reinterpret_cast<char*>(&filename[0]), sizeof(char) * stringSize);
        ifs.read(reinterpret_cast<char*>(&stringSize), sizeof(uint16_t));
        md5->resize(stringSize);
        ifs.read(reinterpret_cast<char*>(&(*md5)[0]), sizeof(char) * stringSize);
        md5s.push_back(std::make_pair(filename,md5));
    }
    
    uLongf decompressedSize;
    ifs.read(reinterpret_cast<char*> (&decompressedSize), sizeof(uint64_t)); // read gz-uncompressed size
    Bytef *data = new Bytef[decompressedSize];
    
    uLong compressedSize = fileSize - sizeof(uint64_t);
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
    
    memcpy(&readLensTmp8[0], ptr, len8 * sizeof(readLensTmp8[0]));
    ptr += len8 * sizeof(uint64_t);
    memcpy(&readLensTmp16[0], ptr, len16 * sizeof(readLensTmp16[0]));
    ptr += len16 * sizeof(uint64_t);
    memcpy(&readLensTmp64[0], ptr, len64 * sizeof(readLensTmp64[0]));
    ptr += len64 * sizeof(uint64_t);
    
    delete[] data;
    
    // add to vector
    readLens.insert(readLensTmp);

    totReads += len8 + len16 + len64;
}

void InReads::printMd5() {
    for (auto md5 : md5s)
        std::cout<<md5.first<<": "<<*md5.second<<std::endl;
}
