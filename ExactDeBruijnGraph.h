#ifndef DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
#define DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H

#include <set>
#include <vector>
#include <unordered_set>
#include "BloomFilter.h"

using namespace std;

class ExactDeBruijnGraph {
public:
    ExactDeBruijnGraph(string inputPath, unsigned int mer_counts, int k);

    void simple_traverse(string inputPath, string outputPath, int maxBreadth, int maxDepth);

private:
    int k;
    BloomFilter bloomFilter;
    unordered_set<string> criticalFP;
    const int M = 100000;

    void initializeBloomFilter(string inputPath);

    void findCriticalFP(string inputPath);

    set<string> findP(set<string> &S);

    bool isPartOfDeBruijnGraph(string kmer);

    string get_lesser(string s);
};

#endif //DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
