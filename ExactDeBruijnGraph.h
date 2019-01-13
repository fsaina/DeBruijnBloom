#ifndef DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
#define DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H

#include <set>
#include <vector>
#include <unordered_set>
#include "BloomFilter.h"

using namespace std;

class ExactDeBruijnGraph {
public:
    ExactDeBruijnGraph(vector<string> &kmers, int k);

    void traverse(vector<string> &kmers, string outputPath, int maxBreadth, int maxDepth);

private:
    int k;
    BloomFilter bloomFilter;
    unordered_set<string> criticalFP;
    const int M = 100000;

    void initializeBloomFilter(vector<string> &kmers);

    void findCriticalFP(vector<string> &kmers);

    vector<string> findP(vector<string> &S);

    bool isPartOfDeBruijnGraph(string kmer);

    string get_lesser(string s);
};

#endif //DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
