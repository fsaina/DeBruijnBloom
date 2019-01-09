#include <iostream>
#include <list>
#include <string.h>
#include <fstream>
#include "ExactDeBruijnGraph.h"
#include "KmerUtil.h"

ExactDeBruijnGraph::ExactDeBruijnGraph(vector<string> &kmers, int k) : k(k), bloomFilter(kmers, k) {
    initializeBloomFilter(kmers);
    findCriticalFP(kmers);
}

// TODO kmers should be loaded from file for the sake of RAM size as described in algorithm
void ExactDeBruijnGraph::initializeBloomFilter(vector<string> &kmers) {
    cout << "Creating Bloom filter..." << endl;
    for (string kmer : kmers) {
        bloomFilter.add(kmer);
//        bloomFilter.add(KmerUtil::reverseComplement(kmer)); TODO add support for reverse complemnts later if have enough RAM
    }
}

// TODO kmers should be loaded from file for the sake of RAM size as described in algorithm, SEQUENTIALLY
// TODO add the reverse complements too if have enough RAM
void ExactDeBruijnGraph::findCriticalFP(vector<string> &kmers) {
    cout << "Finding critical FP set..." << endl;
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
bool ExactDeBruijnGraph::isPartOfDeBruijnGraph(string kmer) {
    return bloomFilter.contains(kmer) && criticalFP.count(kmer) == 0;
}

void ExactDeBruijnGraph::traverse(vector<string> kmers, string outputPath, int maxBreadth, int maxDepth) {
    cout << "Start traversal..." << endl;
    set<string> startingKmers;
    for (string kmer : kmers) {
        vector<string> leftExtensions = KmerUtil::generateLeftExtensions(kmer);

        int inboundCount = 0;
        for (string e : leftExtensions) {
            if (isPartOfDeBruijnGraph(e)) {
                inboundCount++;
            }
        }
        if (inboundCount == 0) {
            startingKmers.insert(kmer);
        }
    }

    unsigned long startingKmersSize = startingKmers.size();
    int kmerIndex = 1; // just for output tracking

    set<string> contigs;
    set<string> marked;
    for (string start : startingKmers) {
        cout << "Starting kmer: " << kmerIndex++ << "/" << startingKmersSize << endl;
        list<string> branches;
        branches.push_back(start);

        int depth = 0;
        while (depth < maxDepth) {
            list<string> branchesToAdd;
            for (string& branch : branches) {
                if (branch.size() - k < depth) continue;

                string lastKmer = KmerUtil::extractLastKmerInSequence(branch, k);
                marked.insert(lastKmer);

                vector<string> extensions = KmerUtil::generateRightExtensions(lastKmer);
                vector<string> validExtensions;
                for (string e : extensions) {
                    if (isPartOfDeBruijnGraph(e) && marked.count(e) == 0) {
                        validExtensions.push_back(e);
                    }
                }
                if (validExtensions.empty()) continue;

                string branchClone(branch);
                unsigned long extensionsSize = validExtensions.size();
                for (int i = 0; i < extensionsSize; ++i) {
                    char &lastChar = validExtensions[i].back();
                    if (i == 0) {
                        branch.push_back(lastChar);
                    } else {
                        string newBranch(branchClone);
                        newBranch.push_back(lastChar);
                        branchesToAdd.push_back(newBranch);
                    }
                }
            }

            for (string branch : branchesToAdd) {
                if (branches.size() >= maxBreadth)
                    break;
                branches.push_back(branch);
            }

            depth++;
        }

        for (string branch : branches) {
            if (branch.length() >= 2*k + 1) {
                contigs.insert(branch);
            }
        }
    }

    ofstream output;
    output.open(outputPath);
    for (string contig : contigs) {
        output << contig << endl;
    }
    output.close();
}
