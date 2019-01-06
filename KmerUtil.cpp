#include "KmerUtil.h"

const string KmerUtil::reverseComplement(string kmer) {
    unsigned long length = kmer.size();

    for (int i = 0; i < length/2; ++i) {
        swap(kmer[i], kmer[length-i-1]);
        kmer[i] = complement(kmer[i]);
        kmer[length-i-1]=complement(kmer[length-i-1]);
    }

    return kmer;
}

char KmerUtil::complement(char c) {
    switch (c) {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        default:
            return 'U';
    }
}
