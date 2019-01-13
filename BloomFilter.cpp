#include "BloomFilter.h"
#include "MurmurHash3.h"
#include <iostream>
#include <string>

using namespace std;

/*
 * Class constructor of a bloom filter
 *
 * Args:
 *      filterSize (int) : Size of the underlying bit array
 *      numHashes (int) : Number of hash functions to use when populating the bit array
 */
BloomFilter::BloomFilter(int filterSize, int numHashes) : bits(filterSize), numHashes(numHashes) {}

/*
 * Class constructor of a bloom filter that calculates the optimal datastructure size depeding
 * on the number of k-mers counts to be stored inside
 *
 * Args:
 *      merCounts (int) : Number of k-mers to be stored in the bloom filter
 *      k (int) : Length of all the mers
 */
BloomFilter::BloomFilter(unsigned int merCounts, int k) {
    int filterSize = calculate_filter_size(merCounts, k);
    int numHashFunctions = calculate_number_of_hash_functions(filterSize, merCounts);

    // adjust filter size
    while (filterSize % numHashFunctions != 0) {
        filterSize++;
    }

    new (this) BloomFilter(filterSize, numHashFunctions);
}

/*
 * Calculates the optimal size of the bit array on the expected number of kmers
 * that will be stored within.
 *
 * Returns optimal filter size
 */
int BloomFilter::calculate_filter_size(int kmerSize, int k) {
    return (int) (1.44 * kmerSize * log2(1/(2.08/(16*k))));
}

/*
 * Calculates the optimal number of hash functions on the expected number of kmers
 * that will be stored within.
 *
 * Returns optimal number of hash functions
 */
int BloomFilter::calculate_number_of_hash_functions(int filterSize, float kmerSize) {
    return (int) filterSize/kmerSize*log(2);
}

/*
 * Map a given string to its appropriate hash number.
 *
 * Returns an array of two hash codes with a size of 64 bits
 */
array<uint64_t, 2> BloomFilter::hash(const string *s) {
    array<uint64_t, 2> hashes;
    uint32_t seed = 420;
    MurmurHash3_x64_128(s->c_str(), (int) s->size(), seed, hashes.data());
    return hashes; // return two 64 bit hashes
}

/*
 * Provided two independent hash values, the total filter size of the bit array and the required n-th
 * hash, return the n-th hash. This implementation of double hashing greatly reduces the number of
 * hash functions required for implementation in this application.
 */
inline uint64_t BloomFilter::nthHash(uint8_t n, uint64_t hashA, uint64_t hashB, uint64_t filterSize) {
    uint64_t h = (hashA + n * hashB) % filterSize;
    return h;
}

/**
 * Add a particular k-mer in the form of a string to the bloom filter
 *
 * @param s - string representation of a k-mer
 */
void BloomFilter::add(const string s) {
    array<uint64_t, 2> hashes = hash(&s);

    for(int n = 0; n < numHashes; n++) {
        bits[nthHash(n, hashes[0], hashes[1], this->bits.size())] = true;
    }
}

/**
 * Test whether of not a particular string is present in the bloom filter.
 * Please note that false positives are possible.
 *
 * @param s - string representation of a k-mer
 * @return boolean if the string is present in the bloom filter.
 */
bool BloomFilter::contains(const string s) {
    array<uint64_t, 2> hashes = hash(&s);

    if (bits.size() == 0) {
        throw "Empty bloom filter array!";
    }

    if (numHashes < 1) {
        throw "Number of hash functions is less than 1!";
    }

    for(int n = 0; n < this->numHashes; n++){
        if (bits[nthHash(n, hashes[0], hashes[1], bits.size())] == false) {
            return false;
        }
    }
    return true;
}
