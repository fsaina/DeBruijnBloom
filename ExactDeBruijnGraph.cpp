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
//        bloomFilter.add(KmerUtil::reverseComplement(kmer)); TODO add support for reverse complemnts later if have enough RAM
    }
}

// TODO kmers should be loaded from file for the sake of RAM size as described in algorithm, SEQUENTIALLY
// TODO add the reverse complements too if have enough RAM
void ExactDeBruijnGraph::findCriticalFP(vector<string> &kmers) {
    set<string> S;
    for (string s : kmers) {
        S.insert(s);
    }

    set<string> P = findP(S);

    for (string p : P) {
        // if S does not contain p
        if (S.find(p) == S.end()) {
            criticalFP.insert(p);
        }
    }
}

set<string> ExactDeBruijnGraph::findP(set<string> &S) {
    set<string> P;

    for (string s : S) {
        vector<string> E = KmerUtil::generateExtensions(s);
        for (string e : E) {
            if (bloomFilter.contains(e)) {
                P.insert(e);
            }
        }

//        TODO reverse complementation
//        vector<string> rcE = KmerUtil::generateExtensions(KmerUtil::reverseComplement(s));
//        for (string e : rcE) {
//            if (bloomFilter.contains(e)) {
//                P.insert(e);
//            }
//        }
    }

    return P;
}

/*
 * Each query to the Bloom filter is modified such that the yes answer is returned
 * if and only if the Bloom filter answers yes and the element is not in cFP.
 */
bool ExactDeBruijnGraph::bloomFilterQuery(string kmer) {
    return bloomFilter.contains(kmer) && criticalFP.count(kmer) == 0;
}

void ExactDeBruijnGraph::traverse(vector<string> kmers, string outputPath, int breadth, int depth) {
    set<string> startingKmers;
    for (string kmer : kmers) {
        vector<string> leftExtensions = KmerUtil::generateLeftExtensions(kmer);

        int inboundCount = 0;
        for (string e : leftExtensions) {
            if (bloomFilterQuery(e)) {
                inboundCount++;
            }
        }
        if (inboundCount == 0) {
            startingKmers.insert(kmer);
        }
    }

    for (string s : startingKmers) {
        cout << s << endl;
    }
    cout << startingKmers.size() << endl;
}
