#include <iostream>
#include "ExactDeBruijnGraph.h"
#include "KmerUtil.h"


ExactDeBruijnGraph::ExactDeBruijnGraph(vector<string> &kmers, int k) : bloomFilter(kmers, k) {
    initializeBloomFilter(kmers);
    findCriticalFP(kmers);
}

// TODO kmers should be loaded from file for the sake of RAM size as described in algorithm
void ExactDeBruijnGraph::initializeBloomFilter(vector<string> &kmers) {
    for (string kmer : kmers) {
        bloomFilter.add(kmer);
        bloomFilter.add(KmerUtil::reverseComplement(kmer));
    }
}

// TODO kmers should be loaded from file for the sake of RAM size as described in algorithm
void ExactDeBruijnGraph::findCriticalFP(vector<string> &kmers) {
    set<string> P = findP(kmers);

    set<string> S;
    for (string s : kmers) {
        S.insert(s);
    }

    for (string p : P) {
        if (S.find(p) == S.end()) { // if S does not contain
            criticalFP.insert(p);
        }
    }
}

set<string> ExactDeBruijnGraph::findP(vector<string> &kmers) {
    // TODO generate extensions(extend suffix and prefix) and check if is not in bloom filter
    return set<string>();
}

/*
 * Each query to the Bloom filter is modified such that the yes answer is returned
 * if and only if the Bloom filter answers yes and the element is not in cFP.
 */
bool ExactDeBruijnGraph::bloomFilterQuery(string kmer) {
    return bloomFilter.contains(kmer) && criticalFP.count(kmer) == 0;
}
