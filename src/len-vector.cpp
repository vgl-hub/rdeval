#include <iterator>
#include <algorithm>
#include <vector>

#include "len-vector.h"


void LenVector::push_back(size_t size) {
    if (size > 255 && size < 65536)
        readLens16.push_back(size);
    else if (size < 256)
        readLens8.push_back(size);
    else
        readLens64.push_back(size);
    
    ++curr;
}

void LenVector::insert(const LenVector &vector) {
    readLens8.insert(std::end(readLens8), std::begin(vector.readLens8), std::end(vector.readLens8));
    readLens16.insert(std::end(readLens16), std::begin(vector.readLens16), std::end(vector.readLens16));
    readLens64.insert(std::end(readLens64), std::begin(vector.readLens64), std::end(vector.readLens64));
    curr += vector.readLens8.size() + vector.readLens16.size() + vector.readLens64.size();
}

uint64_t LenVector::back() {
    
    uint64_t largest = 0;
    
    if (readLens64.size())
        largest = readLens64.back();
    else if (readLens16.size())
        largest = readLens16.back();
    else if (readLens8.size())
        largest = readLens8.back();
    
    return largest;
}

uint64_t LenVector::front() {
    
    uint64_t smallest = 0;
    
    if (readLens8.size())
        smallest = readLens8.front();
    else if (readLens16.size())
        smallest = readLens16.front();
    else if (readLens64.size())
        smallest = readLens64.front();
    
    return smallest;
}

uint64_t LenVector::operator[](uint64_t index) {
    
    uint64_t value;
    
    if (index < 0 || index >= (readLens8.size() + readLens16.size() + readLens64.size()))
        throw std::out_of_range("Index out of bounds");
    
    if (index < readLens8.size())
        value = readLens8[index];
    else if (index < (readLens8.size() + readLens16.size()))
        value = readLens16[index-readLens8.size()];
    else
        value = readLens64[index-readLens8.size()-readLens16.size()];
    
    return value;
}

uint64_t LenVector::size() {
    return (readLens8.size() + readLens16.size() + readLens64.size());
}

void LenVector::sort() {
    std::sort(readLens8.begin(), readLens8.end());
    std::sort(readLens16.begin(), readLens16.end());
    std::sort(readLens64.begin(), readLens64.end());
}

std::vector<uint64_t> LenVector::all() { // far from optimal implementation, but needed without iterator
    std::vector<uint64_t> allReadLens;
    allReadLens.reserve(readLens8.size() + readLens16.size() + readLens8.size()); // preallocate memory
    allReadLens.insert(allReadLens.end(), readLens8.begin(), readLens8.end());
    allReadLens.insert(allReadLens.end(), readLens16.begin(), readLens16.end());
    allReadLens.insert(allReadLens.end(), readLens64.begin(), readLens64.end());
    return allReadLens;
}

std::vector<uint8_t>& LenVector::getReadLens8() {
    return readLens8;
}

std::vector<uint16_t>& LenVector::getReadLens16() {
    return readLens16;
}

std::vector<uint64_t>& LenVector::getReadLens64() {
    return readLens64;
}
