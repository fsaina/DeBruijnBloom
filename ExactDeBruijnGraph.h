#ifndef DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
#define DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H

#include <set>
#include <vector>
#include <unordered_set>
#include "BloomFilter.h"

using namespace std;

/*
 * Our Implementation of ExactDeBruijnGraph
 *
 * Author(s): Marin Vernier, Marin Kukovačec
 */
class ExactDeBruijnGraph {
public:
    ExactDeBruijnGraph(string inputPath, unsigned int mer_counts, int k);
    void simple_traverse(string inputPath, string outputPath, int maxBreadth, int maxDepth);
    void traverse(string inputPath, string outputPath, int maxBreadth, int maxDepth);
    unsigned long graphSizeInBytes();

private:
    int k;
    BloomFilter bloomFilter;
    unordered_set<string> criticalFP;
    const int M = 100000;

    void initializeBloomFilter(string inputPath);
    void findCriticalFP(string inputPath);
    unordered_set<string> findP(unordered_set<string> &S);
    bool isPartOfDeBruijnGraph(string kmer);
    unordered_set<string> loadKmersFromFile(string path);
    bool isSimpleNode(string kmer);
    int countEdgesLeft(string kmer);
    int countEdgesRight(string kmer);
};

#endif //DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
