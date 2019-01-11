#ifndef DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
#define DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H

#include <set>
#include "BloomFilter.h"

using namespace std;

class ExactDeBruijnGraph {
public:
    ExactDeBruijnGraph(string inputPath, unsigned int mer_counts, int k);

    void traverse(string inputPath, string outputPath, int maxBreadth, int maxDepth);

private:
    int k;
    BloomFilter bloomFilter;
    set<string> criticalFP;

    void initializeBloomFilter(string inputPath);

    void findCriticalFP(string inputPath);

    set<string> findP(set<string> &S);

    bool isPartOfDeBruijnGraph(string kmer);
};

#endif //DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
