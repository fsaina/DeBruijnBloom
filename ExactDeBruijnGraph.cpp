#include <iostream>
#include <list>
#include <string.h>
#include <fstream>
#include "ExactDeBruijnGraph.h"
#include "KmerUtil.h"

ExactDeBruijnGraph::ExactDeBruijnGraph(string inputPath, unsigned int mer_counts, int k) : k(k), bloomFilter(mer_counts, k) {
    initializeBloomFilter(inputPath);
    findCriticalFP(inputPath);
}

//  TODO add support for reverse complemnts later if have enough RAM
void ExactDeBruijnGraph::initializeBloomFilter(string inputPath) {
    cout << "Creating Bloom filter..." << endl;

    ifstream in(inputPath);

    string line;
    while(getline(in, line)) {
        if (line.find(">") != 0) {
           bloomFilter.add(line);
        }
    }

    in.close();
}


// TODO add the reverse complements too if have enough RAM
void ExactDeBruijnGraph::findCriticalFP(string inputPath) {
    cout << "Finding critical FP set..." << endl;

    // Create set S, set of all kmers that are in graph
    set<string> S;
    ifstream in(inputPath);

    string line;
    while(getline(in, line)) {
        if (line.find(">") != 0) {
           S.insert(line);
        }
    }

    in.close();

    // Find set P. Set P is set E filtered with Bloom filter. Where E is set of extensions of S(one node extensions).
    set<string> P = findP(S);

    // Set of critical false positives is gained as P \ S.
    for (string p : P) {
        // if S does not contain p
        if (S.find(p) == S.end()) {
            criticalFP.insert(p);
        }
    }
}

/*
 * For each kmer in S, generate extensions(kmers that are neighbours in graph) and add them to P if Bloom filter
 * says that they are part of graph. I.e. set P is set that contains true positives and false positives.
 */
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

void ExactDeBruijnGraph::traverse(string inputPath, string outputPath, int maxBreadth, int maxDepth) {
    cout << "Start traversal..." << endl;

    // Generate a set of starting kmers, i.e. kmers that have zero inbound nodes. They will be starting nodes of contigs.
    set<string> startingKmers;

    ifstream in(inputPath);

    string line;
    while(getline(in, line)) {
        if (line.find(">") != 0) {
          vector<string> leftExtensions = KmerUtil::generateLeftExtensions(line);

          int inboundCount = 0;
          for (string e : leftExtensions) {
              if (isPartOfDeBruijnGraph(e)) {
                  inboundCount++;
              }
          }
          if (inboundCount == 0) {
              startingKmers.insert(line);
          }
        }
    }

    in.close();

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
    for (string contig : contigs) {
        output << contig << endl;
    }
    output.close();
}
