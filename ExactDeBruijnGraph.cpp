#include <iostream>
#include <list>
#include <string.h>
#include <fstream>
#include "ExactDeBruijnGraph.h"
#include "KmerUtil.h"

ExactDeBruijnGraph::ExactDeBruijnGraph(vector<string> &kmers, int k) : k(k), bloomFilter((int) kmers.size(), k) {
    initializeBloomFilter(kmers);
    findCriticalFP(kmers);
}

void ExactDeBruijnGraph::initializeBloomFilter(vector<string> &kmers) {
    cout << "Creating Bloom filter..." << endl;

    for (string kmer : kmers) {
          bloomFilter.add(kmer);
          bloomFilter.add(KmerUtil::reverseComplement(kmer));
    }
}


void ExactDeBruijnGraph::findCriticalFP(vector<string> &kmers) {
    cout << "Finding critical FP set..." << endl;

    // Find set P. Set P is set E filtered with Bloom filter. Where E is set of extensions of S(one node extensions).
    vector<string> P = findP(kmers);

    int k = 0;
    while(k < kmers.size()) {
        unordered_set<string> Pi;
        vector<string> Dn;
        while (Pi.size() < M && k != kmers.size()) {
            Pi.insert(kmers[k]);
            k++;
        }

        for (int i = 0; i != P.size(); ++i) {
            if (Pi.count(get_lesser(P[i])) == 0) {
                Dn.push_back(get_lesser(P[i]));
            }
        }
        P = Dn;
    }
    unordered_set<string> d(P.begin(), P.end());
    criticalFP = d;
}

string ExactDeBruijnGraph::get_lesser(string s) {
    string s2 = KmerUtil::reverseComplement(s);
    return s < s2 ? s : s2;
}

/*
 * For each kmer in S, generate extensions(kmers that are neighbours in graph) and add them to P if Bloom filter
 * says that they are part of graph. I.e. set P is set that contains true positives and false positives.
 */
vector<string> ExactDeBruijnGraph::findP(vector<string> &S) {
    vector<string> P;

    for (string s : S) {
        vector<string> E = KmerUtil::generateExtensions(s);
        for (string e : E) {
            if (bloomFilter.contains(e)) {
                P.push_back(get_lesser(e));
            }
        }
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

/*
 * Function for graph traversal. It starts from a solid k-mer and, using Bloom filter and cFP structure, looks for
 * its neighbours using bounded-breadth, bounded-depth BFS. Method takes the name of the file in which contigs will
 * be saved as an argument.
 */
void ExactDeBruijnGraph::traverse(vector<string> &kmers, string outputPath, int maxBreadth, int maxDepth) {
    cout << "Start traversal..." << endl;

    // Generate a set of starting kmers, i.e. kmers that have zero inbound nodes. They will be starting nodes of contigs.
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

    set<string> contigs; // set of generated contigs
    set<string> marked; // marking structure, used to avoid already seen nodes
    for (string start : startingKmers) {
        cout << "Starting kmer: " << kmerIndex++ << "/" << startingKmersSize << endl;
        list<string> paths;
        paths.push_back(start); // add starting kmer as first path to search

        int depth = 0; // depth is counted in nodes
        while (depth < maxDepth) {
            list<string> pathsToAdd;
            for (string& path : paths) {
                // ignore paths for which search already stopped
                if (path.size() - k < depth) continue;

                string lastKmer = KmerUtil::extractLastKmerInSequence(path, k);
                marked.insert(lastKmer);

                // find neighbours of last kmer in observing path and filter them(check if are part of graph)
                vector<string> extensions = KmerUtil::generateRightExtensions(lastKmer);
                vector<string> validExtensions;
                for (string e : extensions) {
                    if (isPartOfDeBruijnGraph(e) && marked.count(e) == 0) {
                        validExtensions.push_back(e);
                    }
                }

                if (validExtensions.empty()) continue;

                // prolong path with new neighbour nodes
                string pathClone(path);
                unsigned long extensionsSize = validExtensions.size();
                for (int i = 0; i < extensionsSize; ++i) {
                    char &lastChar = validExtensions[i].back();
                    if (i == 0) { // to observing path only add last char of next node
                        path.push_back(lastChar);
                    } else { // if there is more than one valid next node duplicate the current path and branch the search
                        string newPath(pathClone);
                        newPath.push_back(lastChar);
                        pathsToAdd.push_back(newPath);
                    }
                }
            }

            for (string p : pathsToAdd) {
                if (paths.size() >= maxBreadth)
                    break;
                paths.push_back(p);
            }

            depth++;
        }

        for (string p : paths) {
            // contigs shorter that 2k+1 characters are ignored
            if (p.length() >= 2*k + 1) {
                contigs.insert(p);
            }
        }
    }

    ofstream output;
    output.open(outputPath);
    int outputIndex = 0;
    for (string contig : contigs) {
        output << ">" << outputIndex++ << "__len__" << contig.size() << endl;
        output << contig << endl;
    }
    output.close();
}
