#ifndef DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
#define DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H

#include <set>
#include "BloomFilter.h"

using namespace std;

class ExactDeBruijnGraph {
public:
    ExactDeBruijnGraph(vector<string> &kmers, int k);

    void traverse(vector<string> kmers, string outputPath, int maxBreadth, int maxDepth);

private:
    int k;
    BloomFilter bloomFilter;
    set<string> criticalFP;

    void initializeBloomFilter(vector<string> &kmers);

    void findCriticalFP(vector<string> &kmers);

    set<string> findP(set<string> &S);

    bool isPartOfDeBruijnGraph(string kmer);
};

#endif //DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
