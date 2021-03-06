#ifndef DEBRUIJNBLOOM_KMERUTIL_H
#define DEBRUIJNBLOOM_KMERUTIL_H

#include <string>
#include <vector>

using namespace std;

/*
 * Utils definition of functions used in ExtractDeBruijnGraph.cpp
 *
 * Author(s): Marin Vernier, Marin Kukovačec
 */
namespace KmerUtil {
    const vector<string> BASES = { "A", "T", "C", "G" };
    const string reverseComplement(string kmer);
    char complement(char c);
    vector<string> generateExtensions(string kmer);
    vector<string> generateLeftExtensions(string kmer);
    vector<string> generateRightExtensions(string kmer);
    string extractLastKmerInSequence(string sequence, int k);
}

#endif //DEBRUIJNBLOOM_KMERUTIL_H
