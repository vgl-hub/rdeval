#include "reads.h"
#include "input.h"

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
