#ifndef LEN_VECTOR_H
#define LEN_VECTOR_H

#include <iterator>
#include <algorithm>
#include <vector>
#include <numeric>

template<typename T>
class LenVector {

    std::vector<std::pair<uint8_t,T>> readLens8; // up to 255
    std::vector<std::pair<uint16_t,T>> readLens16; // up to 65535
    std::vector<std::pair<uint64_t,T>> readLens64; // max length
    uint64_t curr = 0;
    
public:
    
    void push_back(std::pair<uint64_t,T> obj);
    
    void insert(const LenVector &vector);
    
    uint64_t back();
    
    uint64_t front();
    
    std::pair<uint64_t,T> operator[](uint64_t index);
    
    uint64_t size();
    
    void sort();
    
    std::vector<std::pair<uint64_t,T>> all();
    
    std::vector<std::pair<uint8_t,T>>& getReadLens8();

    std::vector<std::pair<uint16_t,T>>& getReadLens16();
    
    std::vector<std::pair<uint64_t,T>>& getReadLens64();
    
    std::vector<uint64_t> getSortedIndex();
};

template<typename T>
void LenVector<T>::push_back(std::pair<uint64_t,T> obj) {
    size_t size = obj.first;
    if (size > 255 && size < 65536) {
        std::pair<uint16_t, T> readObj(size, obj.second);
        readLens16.push_back(readObj);
    }else if (size < 256) {
        std::pair<uint8_t, T> readObj(size, obj.second);
        readLens8.push_back(readObj);
    }else {
        std::pair<uint64_t, T> readObj(size, obj.second);
        readLens64.push_back(readObj);
    }
    ++curr;
}

template<typename T>
void LenVector<T>::insert(const LenVector &vector) {
    readLens8.insert(std::end(readLens8), std::begin(vector.readLens8), std::end(vector.readLens8));
    readLens16.insert(std::end(readLens16), std::begin(vector.readLens16), std::end(vector.readLens16));
    readLens64.insert(std::end(readLens64), std::begin(vector.readLens64), std::end(vector.readLens64));
    curr += vector.readLens8.size() + vector.readLens16.size() + vector.readLens64.size();
}

template<typename T>
uint64_t LenVector<T>::back() {
    
    uint64_t smallest = 0;
    
    if (readLens8.size())
        smallest = readLens8.back().first;
    else if (readLens16.size())
        smallest = readLens16.back().first;
    else if (readLens64.size())
        smallest = readLens64.back().first;
    
    return smallest;
}

template<typename T>
uint64_t LenVector<T>::front() {
    
    uint64_t largest = 0;
    
    if (readLens64.size())
        largest = readLens64.front().first;
    else if (readLens16.size())
        largest = readLens16.front().first;
    else if (readLens8.size())
        largest = readLens8.front().first;
    
    return largest;
}

template<typename T>
std::pair<uint64_t,T> LenVector<T>::operator[](uint64_t index) {
    
    if (index >= (readLens8.size() + readLens16.size() + readLens64.size()))
        throw std::out_of_range("Index out of bounds");
    
    if (index < readLens8.size())
        return readLens8[index];
    else if (index < readLens16.size())
        return readLens16[index-readLens8.size()];
    else
        return readLens64[index-readLens8.size()-readLens16.size()];
}

template<typename T>
uint64_t LenVector<T>::size() {
    return (readLens8.size() + readLens16.size() + readLens64.size());
}

template <class T1, class T2>
struct sort_pair_first {
    bool operator()(const std::pair<T1,T2>&a, const std::pair<T1,T2>&b) {
        if (a.first != b.first)
            return a.first > b.first; // Sort by first element descending
        else
            return a.second > b.second; // If first elements are equal, sort by second element descending
    }
};

template<typename T>
void LenVector<T>::sort() { // sort reads by length ([0])
    
    std::sort(readLens8.begin(), readLens8.end(), sort_pair_first<uint8_t, T>());
    std::sort(readLens16.begin(), readLens16.end(), sort_pair_first<uint16_t, T>());
    std::sort(readLens64.begin(), readLens64.end(), sort_pair_first<uint64_t, T>());
}

template<typename T>
std::vector<std::pair<uint64_t,T>> LenVector<T>::all() { // far from optimal implementation, but needed without iterator
    std::vector<std::pair<uint64_t,T>> allReadLens;
    allReadLens.reserve(readLens8.size() + readLens16.size() + readLens8.size()); // preallocate memory
    allReadLens.insert(allReadLens.end(), readLens8.begin(), readLens8.end());
    allReadLens.insert(allReadLens.end(), readLens16.begin(), readLens16.end());
    allReadLens.insert(allReadLens.end(), readLens64.begin(), readLens64.end());
    return allReadLens;
}

template<typename T>
std::vector<std::pair<uint8_t,T>>& LenVector<T>::getReadLens8() {
    return readLens8;
}

template<typename T>
std::vector<std::pair<uint16_t,T>>& LenVector<T>::getReadLens16() {
    return readLens16;
}

template<typename T>
std::vector<std::pair<uint64_t,T>>& LenVector<T>::getReadLens64() {
    return readLens64;
}

template<typename T>
std::vector<uint64_t> LenVector<T>::getSortedIndex() {
    
    // sort note: range and zip c++20 extension would work well here
    std::vector<uint64_t> idx(size());
    iota(idx.begin(), idx.end(), 0);
    
    // sort indexes based on comparing values in v using std::stable_sort instead of std::sort to avoid unnecessary index re-orderings when v contains elements of equal values
    stable_sort(idx.begin(), idx.begin()+readLens8.size(), [this](size_t i1, size_t i2) { return readLens8[i1] < readLens8[i2]; });
    stable_sort(idx.begin()+readLens8.size(), idx.begin()+readLens16.size(), [this](size_t i1, size_t i2) { return readLens16[i1] < readLens16[i2]; });
    stable_sort(idx.begin()+readLens16.size(), idx.begin()+readLens64.size(), [this](size_t i1, size_t i2) { return readLens64[i1] < readLens64[i2]; });
    
    return idx;
}

#endif /* LEN_VECTOR_H */
