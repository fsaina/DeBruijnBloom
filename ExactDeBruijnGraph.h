#ifndef DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
#define DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H

#include <set>
#include "BloomFilter.h"

using namespace std;

class ExactDeBruijnGraph {
public:
    ExactDeBruijnGraph(vector<string> &kmers, int k);

    // TODO graph traversal function

private:
    BloomFilter bloomFilter;
    set<string> criticalFP;

    void initializeBloomFilter(vector<string> &kmers);

    void findCriticalFP(vector<string> &kmers);

    set<string> findP(set<string> &S);

    bool bloomFilterQuery(string kmer);
};

#endif //DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
