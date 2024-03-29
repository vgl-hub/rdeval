#include <vector>
#include <string>
#include <fstream>
#include <memory>
#include <condition_variable>
#include <mutex>

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

void Input::read(InReads& inReads) {
    
    if (userInput.inReads.empty()) {return;}

    inReads.load(&userInput);

}
