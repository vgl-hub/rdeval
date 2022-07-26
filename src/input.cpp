#include <stdlib.h>
#include <vector>
#include <string>
#include <queue>
#include <thread>

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
#include "threadpool.h"
#include "reads.h"
#include "stream-obj.h"

#include "input.h"

void Input::load(UserInput userInput) {
    
    this->userInput = userInput;
    
}

void Input::read(InReads& inReads) {
    
    if (userInput.iReadFileArg.empty()) {return;}
        
    threadPool.init(maxThreads); // initialize threadpool

    inReads.load(userInput);
    
    while (true) {
        
        if (threadPool.empty()) {threadPool.join(); break;}
        lg.verbose("Remaining jobs: " + std::to_string(threadPool.queueSize()), true);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        
    }

}
