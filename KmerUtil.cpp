#include "KmerUtil.h"

const string KmerUtil::reverseComplement(string kmer) {
    unsigned long length = kmer.size();

    for (int i = 0; i <= length/2; ++i) {
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

vector<string> KmerUtil::generateExtensions(string kmer) {
    vector<string> E;

    unsigned long kmerSize = kmer.size();
    for (string b : BASES) {
        // left extensions
        E.push_back(b + kmer.substr(0, kmerSize - 1));
        // right extension
        E.push_back(kmer.substr(1) + b);
    }

    return E;
}
