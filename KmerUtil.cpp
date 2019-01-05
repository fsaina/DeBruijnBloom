#include "KmerUtil.h"

const string KmerUtil::reverseComplement(string kmer) {
    unsigned long length = kmer.size();
    string reversed(length, '-');

    for (int i = 0; i < length; ++i) {
        reversed[i] = complement(kmer[length - i - 1]);
    }

    return reversed;
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