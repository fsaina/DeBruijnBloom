#include "measures.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

/**
 * Returns the N50 value of the passed list of numbers.
 *
 * Based on the Broad Institute definition:
 * https://www.broad.harvard.edu/crd/wiki/index.php/N50
 *
 * @param lengths - vector of sequence lengths
 * @return N50 measure value for the input lengths
 *
 * Author(s): Marin Vernier
 */
float measures::n50(vector<int> lengths) {
    map<int, int> freq;

    for (int i : lengths) {
        freq[i] = freq[i] + 1;
    }

    vector<int> tmpLengths;

    for (const auto& p : freq) {
        int key = p.first;
        int value = p.second;

        for (int i=0; i<value*key; i++) {
            tmpLengths.push_back(key);
        }
    }

    size_t size = tmpLengths.size();

    if (size % 2 == 0) {
        return (tmpLengths[size / 2 - 1] + tmpLengths[size / 2]) / 2.;
    } else {
        return tmpLengths[size / 2];
    }
}

/*
 * Read a fasta file and retrieve the sequence counts (lengths) assigned to
 * each sequence
 */
vector<int> measures::read_seq_counts(string inputPath){
    ifstream in(inputPath);

    vector<int> counts;
    string delimiter = "__len__";

    string line;
    while(getline(in, line)) {
        if (line.find(">") == 0) {
            string token = line.substr(line.find(delimiter) + delimiter.size());
            counts.push_back(stoi(token));
        }
    }

    in.close();
    return counts;
}
