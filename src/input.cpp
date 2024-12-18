#include <vector>
#include <string>
#include <fstream>
#include <memory>
#include <condition_variable>
#include <mutex>
#include <deque>

#include "global.h"
#include "log.h"
#include "bed.h"
#include "struct.h"
#include "gfa-lines.h"
#include "uid-generator.h"
#include "gfa.h"
#include "functions.h" // global functions
#include "reads.h"
#include "stream-obj.h"

#include "input.h"

void Input::load(UserInputRdeval userInput) {
    
    this->userInput = userInput;
    
}

void Input::read() {
    
    if (userInput.inFiles.empty()) {return;}

    std::string outFile = "";
    if(userInput.outFiles.size())
        outFile = userInput.outFiles[0];
    
    InReads inReads(userInput, outFile); // initialize sequence collection object
    lg.verbose("Read object generated");
    threadPool.init(maxThreads); // initialize threadpool
    
    if (inReads.isOutputBam())
        inReads.writeBamHeader();
    
    inReads.load();
    jobWait(threadPool);
    inReads.writeToStream(); // write last batch
    
    lg.verbose("Data loaded");
    
    if (inReads.isOutputBam())
        inReads.closeBam();
    
    lg.verbose("Generating output");
    
    if (userInput.md5_flag)
       inReads.printMd5();
    else if (userInput.stats_flag) // output summary statistics
        inReads.report();
    else if (userInput.outSize_flag)
        inReads.printReadLengths();
    else if (userInput.quality_flag)
        inReads.printQualities();
    else if (userInput.content_flag)
        inReads.printContent();
    
    if (userInput.outFiles.size()) {
        lg.verbose("Writing rd file(s)");
        for (std::string file : userInput.outFiles)
            if (getFileExt(file) == "rd")
                inReads.printTableCompressed(file);
    }
}
