#ifndef DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
#define DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H

#include <set>
#include "BloomFilter.h"

using namespace std;

class ExactDeBruijnGraph {
public:
    ExactDeBruijnGraph(vector<string> &kmers, int k);

private:
    BloomFilter bloomFilter;
    set<string> criticalFP;

    void findCriticalFP(vector<string> &kmers);
};

#endif //DEBRUIJNBLOOM_EXACTDEBRUIJNGRAPH_H
