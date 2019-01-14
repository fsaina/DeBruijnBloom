#include <iostream>
#include <list>
#include <string.h>
#include <fstream>
#include "ExactDeBruijnGraph.h"
#include "KmerUtil.h"

/*
 * Class constructor that initializes bloom filter and finds critical false positives
 *
 * Author(s): Marin Vernier, Marin Kukovačec
 */
ExactDeBruijnGraph::ExactDeBruijnGraph(string inputPath, unsigned int mer_counts, int k) : k(k), bloomFilter(mer_counts, k) {
    initializeBloomFilter(inputPath);
    findCriticalFP(inputPath);
}

/*
 * Loads kmers from file given in arguments and populates the Bloom filter with them.
 *
 * Author(s): Marin Vernier, Marin Kukovačec
 */
void ExactDeBruijnGraph::initializeBloomFilter(string inputPath) {
    cout << "Creating Bloom filter..." << endl;

    unordered_set<string> kmers = loadKmersFromFile(inputPath);

    for (string kmer : kmers) {
        bloomFilter.add(kmer);
    }
}

/*
 * Function for computation of critical false positives which is formally defined as cFP = P\S
 *
 * Author(s): Marin Vernier, Marin Kukovačec
 */
void ExactDeBruijnGraph::findCriticalFP(string inputPath) {
    cout << "Finding critical FP set..." << endl;

    // Create set S, set of all kmers that are in graph
    // TODO try to make sequential read from algorithm in paper, but dont use reverse complements
    unordered_set<string> S = loadKmersFromFile(inputPath);

    // Find set P. Set P is set E filtered with Bloom filter. Where E is set of extensions of S(one node extensions).
    unordered_set<string> P = findP(S);

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
 *
 * Author(s): Marin Vernier, Marin Kukovačec
 */
unordered_set<string> ExactDeBruijnGraph::findP(unordered_set<string> &S) {
    unordered_set<string> P;

    for (string s : S) {
        vector<string> E = KmerUtil::generateExtensions(s);
        for (string e : E) {
            if (bloomFilter.contains(e)) {
                P.insert(e);
            }
        }
    }

    return P;
}

/*
 * Each query to the Bloom filter is modified such that the yes answer is returned
 * if and only if the Bloom filter answers yes and the element is not in cFP.
 *
 * Author(s): Marin Vernier
 */
bool ExactDeBruijnGraph::isPartOfDeBruijnGraph(string kmer) {
    return bloomFilter.contains(kmer) && criticalFP.count(kmer) == 0;
}

/*
 * Method used to print graph size. Size of Bloom filter + size of cFP structure.
 *
 * Author(s): Marin Vernier
 */
unsigned long ExactDeBruijnGraph::graphSizeInBytes() {
    unsigned long size = 0;
    size += bloomFilter.sizeInBytes();
    size += criticalFP.size() * k;
    return size;
}

/*
 * Function for graph traversal. It starts from a solid k-mer and, using Bloom filter and cFP structure, looks for
 * its neighbours using bounded-breadth, bounded-depth BFS. Method takes the name of the file in which contigs will
 * be saved as an argument.
 *
 * Author(s): Marin Vernier
 */
void ExactDeBruijnGraph::simple_traverse(string inputPath, string outputPath, int maxBreadth, int maxDepth) {
    cout << "Start traversal..." << endl;

    // Generate a set of starting kmers, i.e. kmers that have zero inbound nodes. They will be starting nodes of contigs.
    unordered_set<string> startingKmers;
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

    unordered_set<string> contigs; // set of generated contigs
    unordered_set<string> marked; // marking structure, used to avoid already seen nodes
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

/*
 * Function for graph traversal. It starts from a solid k-mer and, using Bloom filter and cFP structure, looks for
 * its neighbours using bounded-breadth, bounded-depth BFS. Method takes the name of the file in which contigs will
 * be saved as an argument.
 *
 * Author(s): Marin Vernier
 */
void ExactDeBruijnGraph::traverse(string inputPath, string outputPath, int maxBreadth, int maxDepth) {
    cout << "Start traversal..." << endl;

    // Generate a set of starting kmers, i.e. kmers that have zero inbound nodes. They will be starting nodes of contigs.
    unordered_set<string> startingKmers;
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

    unordered_set<string> contigs; // set of generated contigs
    unordered_set<string> marked; // marking structure, used to avoid already seen nodes
    for (string start : startingKmers) {
        cout << "Starting kmer: " << kmerIndex++ << "/" << startingKmersSize << endl;
        list<string> paths;
        paths.push_back(start); // add starting kmer as first path to search

        int depth = 0; // depth is counted in nodes
        while (paths.size() > 0) {
            if (depth >= maxDepth) { // if BFS search reaches the max value bound
                contigs.insert(paths.back());
                break;
            }

            // auxiliary lists used to memorize elements, for the sake of not updating the iterating list
            list<string> pathsToAdd;
            list<string> pathsToRemove;

            int breadth = paths.size();
            // if path is simple(no branching) just expand the current contig
            if (breadth == 1) {
                for (string& path : paths) {
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

                    unsigned long validExtensionsCount = validExtensions.size();
                    if (validExtensionsCount == 0) {
                        // end of search
                        contigs.insert(path);
                        pathsToRemove.push_back(path);
                    } else if (validExtensionsCount == 1) {
                        // expand current contig
                        char &lastChar = validExtensions[0].back();
                        path.push_back(lastChar);
                    } else {
                        // here starts the branch, save current contig and make branches, from now BFS will start
                        contigs.insert(path);
                        pathsToRemove.push_back(path);
                        for(string validE : validExtensions) {
                            pathsToAdd.push_back(validE);
                        }
                        depth++;
                    }
                }
            } else {
                // there are more active branches, BFS search
                depth++;
                for (string& path : paths) {
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

                    string pathClone(path);
                    unsigned long validExtensionsCount = validExtensions.size();
                    if (validExtensionsCount == 0) {
                        // current branch is over, this way if one branch remains active the search is going to continue normally
                        pathsToRemove.push_back(path);
                    } else if (validExtensionsCount == 1) {
                        // expand current branch
                        char &lastChar = validExtensions[0].back();
                        path.push_back(lastChar);
                    } else {
                        // try to branch more
                        for(string validE : validExtensions) {
                            char &lastChar = validE.back();
                            string newPath(pathClone);
                            newPath.push_back(lastChar);
                            pathsToAdd.push_back(newPath);
                        }
                    }
                }
            }

            for (string p : pathsToRemove) {
                paths.remove(p);
            }

            if (paths.size() == 1) {
                depth = 0;
            }

            for (string p : pathsToAdd) {
                if (paths.size() >= maxBreadth)
                    break;
                paths.push_back(p);
            }
        }
    }

    ofstream output;
    output.open(outputPath);
    int outputIndex = 0;
    for (string contig : contigs) {
        // contigs shorter that 2k+1 characters are ignored
        if (contig.length() < 2*k + 1) {
            continue;
        }
        output << ">" << outputIndex++ << "__len__" << contig.size() << endl;
        output << contig << endl;
    }
    output.close();
}

/*
 * Loads kmers from file into and unordered set.
 *
 * Author(s): Marin Kukovačec
 */
unordered_set<string> ExactDeBruijnGraph::loadKmersFromFile(string path) {
    ifstream in(path);

    unordered_set<string> S;
    string line;
    while(getline(in, line)) {
        if (line.find(">") != 0) {
            S.insert(line);
        }
    }

    in.close();

    return S;
}

/*
 * Checks if given kmer is simple, i.e. has at most one inbound and one outbound edge.
 *
 * Author(s): Marin Vernier
 */
bool ExactDeBruijnGraph::isSimpleNode(string kmer) {
    if (countEdgesLeft(kmer) > 1)
        return false;
    if (countEdgesRight(kmer) > 1)
        return false;
    return true;
}

/*
 * Count number of possible inbound edges.
 *
 * Author(s): Marin Vernier
 */
int ExactDeBruijnGraph::countEdgesLeft(string kmer) {
    vector<string> extensions = KmerUtil::generateLeftExtensions(kmer);
    int count = 0;
    for (string extension : extensions) {
        if (isPartOfDeBruijnGraph(extension))
            count++;
    }
    return count;
}

/*
 * Count number of possible outbound edges.
 *
 * Author(s): Marin Vernier
 */
int ExactDeBruijnGraph::countEdgesRight(string kmer) {
    vector<string> extensions = KmerUtil::generateRightExtensions(kmer);
    int count = 0;
    for (string extension : extensions) {
        if (isPartOfDeBruijnGraph(extension))
            count++;
    }
    return count;
}
