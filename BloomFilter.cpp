#include "BloomFilter.h"
#include "MurmurHash3.h"
#include <iostream>
#include <string>

using namespace std;

BloomFilter::BloomFilter(int filterSize, int numHashes) : numHashes(numHashes), bits(filterSize) {}

BloomFilter::BloomFilter(vector<string>& kmers, int k) {
    int filterSize = calculate_filter_size((int) kmers.size(), k);
    int numHashFunctions = calculate_number_of_hash_functions(filterSize, (int) kmers.size());

    // adjust filter size
    while (filterSize % numHashFunctions != 0) {
        filterSize++;
    }

    new (this) BloomFilter(filterSize, numHashFunctions);

    for (string s : kmers) {
        this->add(s);
    }
}

int BloomFilter::calculate_filter_size(int kmerSize, int k) {
    return (int) (1.44 * kmerSize * log2(1/(2.08/(16*k))));
}

int BloomFilter::calculate_number_of_hash_functions(int filterSize, float kmerSize) {
    return (int) filterSize/kmerSize*log(2);
}

std::array<uint64_t, 2> BloomFilter::hash(const std::string *s) {
    array<uint64_t, 2> hashes;
    uint32_t seed = 420;
    MurmurHash3_x64_128(s, (int) s->size(), seed, hashes.data());
    return hashes; // return two 64 bit hashes
}

inline uint64_t BloomFilter::nthHash(uint8_t n, uint64_t hashA, uint64_t hashB, uint64_t filterSize) {
    // https://en.wikipedia.org/wiki/Double_hashing
    uint64_t h = (hashA + n * hashB) % filterSize;
    return h;
}

void BloomFilter::add(const string & s) {
    array<uint64_t, 2> hashes = hash(&s);

    for(int n = 0; n < this->numHashes; n++) {
        this->bits[nthHash(n, hashes[0], hashes[1], this->bits.size())] = true;
    }
}

bool BloomFilter::contains(const string & s) {
    array<uint64_t, 2> hashes = hash(&s);

    if (this->bits.size() == 0) {
        throw "Empty bloom filter array!";
    }

    if (numHashes < 1) {
        throw "Number of hash functions is less than 1!";
    }

    for(int n = 0; n < this->numHashes; n++){
        if (this->bits[nthHash(n, hashes[0], hashes[1], bits.size())] == false) {
            return false;
        }
    }
    return true;
}
