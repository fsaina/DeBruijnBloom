#include "ExactDeBruijnGraph.h"

ExactDeBruijnGraph::ExactDeBruijnGraph(vector<string> &kmers, int k) : bloomFilter(kmers, k) {
    findCriticalFP(kmers);
}

void ExactDeBruijnGraph::findCriticalFP(vector<string> &kmers) {

}

