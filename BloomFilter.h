//
// Custom implementation followed http://blog.michaelschmatz.com/2016/04/11/how-to-write-a-bloom-filter-cpp/
//

#ifndef DEBRUJINBLOOM_BLOOMFILTER_H
#define DEBRUJINBLOOM_BLOOMFILTER_H

#include "vector"
#include "array"
#include <math.h>

using namespace std;

class BloomFilter {
    int numHashes;
    vector<bool> bits;

public:
    BloomFilter(int size, int numHashes);
    BloomFilter(vector<string> &kmers, int k);
    void add(const string & s);
    bool contains(const string & s);

private:
    array<uint64_t, 2> hash(const std::string *s);
    uint64_t nthHash(uint8_t n, uint64_t hashA, uint64_t hashB, uint64_t filterSize);
    int calculate_filter_size(int kmerSize, int k);
    int calculate_number_of_hash_functions(int filterSize, float kmerSize);
};

#endif //DEBRUJINBLOOM_BLOOMFILTER_H
