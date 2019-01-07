#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

#include "cmdline.h"
#include "Tests.h"

using namespace std;

static const int K_MER_THRESHOLD = 4;

vector<string> read_file_in_vector(string inputPath){
    vector<string> mers;

    ifstream in(inputPath);

    string line;
    while(getline(in, line)) {
       if (line.size() > K_MER_THRESHOLD) {
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
    p.add<string>("tmp", 't', "Path to directory where to write temporary files", false, "./tmp");

    p.parse_check(argc, argv);

    int k = p.get<int>("kmers");
    int minAbundance = p.get<int>("minAbundance");
    string inputPath = p.get<string>("input");
    string outputPath = p.get<string>("output");
    string jellyfishBinPath = p.get<string>("jellyfish");
    string tmpDir = p.get<string>("tmp");

    string jellyfishTmpFileName = "tmp.fa";
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
    command = jellyfishBinPath + " dump " + tmpDir + defaultJellyfishOutput + "_0 > " + tmpDir + jellyfishTmpFileName;
    system(command.c_str());

    // read the file and load only the k-mers (not their counts)
    vector<string> kmers = read_file_in_vector(tmpDir + jellyfishTmpFileName);

    Tests::run_all_tests(kmers, k);

    return 0;
}
