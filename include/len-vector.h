#ifndef LEN_VECTOR_H
#define LEN_VECTOR_H


class LenVector {

    std::vector<uint8_t> readLens8; // up to 255
    std::vector<uint16_t> readLens16; // up to 65535
    std::vector<uint64_t> readLens64; // max length
    uint64_t curr = 0;
    
public:
    
    void push_back(size_t size);
    
    void insert(const LenVector &vector);
    
    uint64_t back();
    
    uint64_t front();
    
    uint64_t operator[](uint64_t index);
    
    uint64_t size();
    
    void sort();
    
    std::vector<uint64_t> all();
    
    std::vector<uint8_t>& getReadLens8();

    std::vector<uint16_t>& getReadLens16();
    
    std::vector<uint64_t>& getReadLens64();
};

#endif /* LEN_VECTOR_H */
