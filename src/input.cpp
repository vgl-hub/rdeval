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
#include "reads.h"
#include "stream-obj.h"

#include "input.h"

void Input::load(UserInputRdeval userInput) {
    
    this->userInput = userInput;
    
}

void Input::read() {
    
    if (userInput.inFiles.empty()) {return;}

    InReads inReads(userInput); // initialize sequence collection object
    lg.verbose("Read object generated");
    threadPool.init(maxThreads); // initialize threadpool
    inReads.load();
    jobWait(threadPool);

    if (userInput.stats_flag) // output summary statistics
        inReads.report();
    else if (userInput.outSize_flag)
        inReads.printReadLengths();
    else if (userInput.quality_flag)
        inReads.printQualities();
    else if (userInput.content_flag)
        inReads.printContent();
}
