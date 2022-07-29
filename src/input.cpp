#include <stdlib.h>
#include <vector>
#include <string>
#include <queue>
#include <thread>
#include <cstring>

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

    inReads.load(userInput);

}
