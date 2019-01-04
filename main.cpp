#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

//#include "BloomFilter.h"
#include "StringBloomFilter.h"

using namespace std;

vector<string> read_file_in_vector(string inputPath){
    vector<string> mers;

    ifstream in(inputPath);

    string line;
    while(getline(in, line)) {
       if (line.size() > 4) {  // TODO fix magic number
           mers.push_back(line);
       }
    }

    in.close();
    return mers;
}

int main() {
    int k = 21;
    int minAbundance = 3;

    string inputPath = "/Users/filipsaina/CLionProjects/DeBrujinBloom/data/ecoli.fasta";
    string outputPath = "";

    string workingDir = "/Users/filipsaina/CLionProjects/DeBrujinBloom/";
    string tmpDir = workingDir + "tmp/";

    string jellyfishBinPath = workingDir + "jellyfish";
    string jellyfishTmpFilePath = "tmp.fa";
    string defaultJellyfishOutput = "mer_counts.jf";

    string command = "";

    // Remove and create tmp directory
    command = "rm -rf " + tmpDir;
    system(command.c_str());

    command = "mkdir " + tmpDir;
    system(command.c_str());

    // Count all k-mers https://github.com/gmarcais/Jellyfish/tree/master/doc
    command = jellyfishBinPath + " count -m " + to_string(k)
              + " -s 100M -t 4 -C -L " + to_string(minAbundance)
              + " " + inputPath
              + " -o " + tmpDir + defaultJellyfishOutput;
    system(command.c_str());

    // Convert binary jellyfish output to human readable format
    command = jellyfishBinPath + " dump " + tmpDir + defaultJellyfishOutput + "_0" " > " + tmpDir + jellyfishTmpFilePath;
    system(command.c_str());

    // read the file and load only the k-mers (not their counts)
    vector<string> kmers = read_file_in_vector(tmpDir + jellyfishTmpFilePath);

    BloomFilter bf = BloomFilter(kmers, k);


    // TODO remove these tests (all code until return)

    // NOTE test that there are no false negatives
    for (string s : kmers) {
        bool contains = bf.contains(s);
        if (contains == false) {
            printf("> ERROR\n");
        }
    }

    bool contains;
    // NOTE in bloom filter data structure false positives are possible, but false negatives should NEVER happen

    contains = bf.contains("AAAAAAAAAAAAAAAAAAAAA");
    cout << "Should be true : " << boolalpha  << contains << endl;

    contains = bf.contains("ATGAACGGCAGCGGGCCAAAA");
    cout << "Should be false : " << boolalpha << contains << endl;

    return 0;
}
