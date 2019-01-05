#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

#include "cmdline.h"
#include "BloomFilter.h"

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

int main(int argc, char *argv[]) {
    cmdline::parser p;

    p.add<int>("kmers", 'k', "Number of k mers", true);
    p.add<int>("minAbundance", 'a', "Minimal abundance of sequences required", true);
    p.add<string>("input", 'i', "Path to the input fasta file", true);
    p.add<string>("output", 'o', "Path to directory where to write the results of execution", false, "./");
    p.add<string>("jellyfish", 'j', "Path to jellyfish executable", false, "./jellyfish");

    p.parse_check(argc, argv);

    int k = p.get<int>("kmers");
    int minAbundance = p.get<int>("minAbundance");
    string inputPath = p.get<string>("input");
    string outputPath = p.get<string>("output");
    string jellyfishBinPath = p.get<string>("jellyfish");

    string workingDir = "./";
    string tmpDir = workingDir + "tmp/";

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
