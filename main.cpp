#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include "BloomFilter.h"

using namespace std;

void bloomDemo();

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

    string inputPath = "./data/ecoli.fasta";
    string outputPath = "";

    string workingDir = "./";
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
    command = jellyfishBinPath + " dump " + tmpDir + defaultJellyfishOutput + " > " + tmpDir + jellyfishTmpFilePath;
    system(command.c_str());

    // read the file and load only the k-mers (not their counts)
    vector<string> mers = read_file_in_vector(tmpDir + jellyfishTmpFilePath);

    // TODO remove me later
    for (string s : mers) {
        cout << s << endl;
    }

    // TODO remove this function call
    bloomDemo();

    return 0;
}

/*
 * TODO remove this function
 */
void bloomDemo() {
    BloomFilter bloomFilter;

    bloomFilter.add("kuki");
    bloomFilter.add("sajo");
    bloomFilter.add("verni");

    cout << "Should be true: " << (bloomFilter.exists("kuki") ? "true" : "false") << endl;
    cout << "Should be false: " << (bloomFilter.exists("brek") ? "true" : "false") << endl;
}
