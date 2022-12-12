#ifndef INPUT_H
#define INPUT_H

class Input {
    
    UserInputRdeval userInput;
    
    //intermediates
    std::string h;
    
    // stream read variable definition
    std::string firstLine;
    
    StreamObj streamObj;
    
    std::string newLine, seqHeader, seqComment, line;
    
    std::shared_ptr<std::istream> stream;
    
public:
    
    void load(UserInputRdeval userInput);
    
    void read(InReads& inReads);
    
};

#endif /* INPUT_H */
