#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

#include "cmdline.h"
#include "Tests.h"
#include "ExactDeBruijnGraph.h"

using namespace std;

vector<string> read_file_in_vector(string inputPath){
    vector<string> mers;

    ifstream in(inputPath);

    string line;
    while(getline(in, line)) {
        // ignore lines starting with '>'
        if (line.find(">") != 0) {
           mers.push_back(line);
        }
    }

    in.close();
    return mers;
}

int main(int argc, char *argv[]) {
    clock_t start = clock();

    cmdline::parser p;

    p.add<int>("kmers", 'k', "Number of k mers", true);
    p.add<int>("minAbundance", 'a', "Minimal abundance of sequences required", true);
    p.add<int>("maxBreadth", 'b', "Maximum node breadth while doing a BFS traversal of the graph", true);
    p.add<int>("maxDepth", 'd', "Maximum node breadth while doing a BFS traversal of the graph", true);
    p.add<string>("input", 'i', "Path to the input fasta file", true);
    p.add<string>("output", 'o', "Path to directory where to write the results of execution", false, "./output");
    p.add<string>("jellyfish", 'j', "Path to jellyfish executable", false, "./bin/jellyfish");
    p.add<string>("tmp", 't', "Path to directory where to write temporary files", false, "./tmp");

    p.parse_check(argc, argv);

    int k = p.get<int>("kmers");
    int minAbundance = p.get<int>("minAbundance");
    int maxBreadth = p.get<int>("maxBreadth");
    int maxDepth = p.get<int>("maxDepth");
    string inputPath = p.get<string>("input");
    string outputPath = p.get<string>("output");
    string jellyfishBinPath = p.get<string>("jellyfish");
    string tmpDir = p.get<string>("tmp");

    string jellyfishTmpFileName = "tmp.fa";
    string defaultJellyfishOutput = "mer_counts.jf";
    string defaultProgramOutput = "output.fasta";

    string command = "";

    // Remove and create tmp directory
    command = "rm -rf " + tmpDir;
    system(command.c_str());

    command = "mkdir " + tmpDir;
    system(command.c_str());

    // Count all k-mers https://github.com/gmarcais/Jellyfish/tree/master/doc
    command = jellyfishBinPath + " count -m " + to_string(k)
              + " -s 100M -t 4 -L " + to_string(minAbundance) // TODO maybe return '-C' flag
              + " " + inputPath
              + " -o " + tmpDir + '/' + defaultJellyfishOutput;
    system(command.c_str());

    // Convert binary jellyfish output to human readable format
    command = jellyfishBinPath + " dump " + tmpDir + '/' + defaultJellyfishOutput + " > " + tmpDir + '/' + jellyfishTmpFileName;
    system(command.c_str());

    // read the file and load only the k-mers (not their counts)
    vector<string> kmers = read_file_in_vector(tmpDir + '/' + jellyfishTmpFileName);

    cout<< "Number of kmers: " << kmers.size()<< endl;

    command = "rm -rf " + outputPath;
    system(command.c_str());

    command = "mkdir " + outputPath;
    system(command.c_str());

    ExactDeBruijnGraph graph = ExactDeBruijnGraph(kmers, k);
    graph.traverse(kmers, outputPath + "/" + defaultProgramOutput, maxBreadth, maxDepth);

//    Tests::run_all_tests(kmers, k);

    clock_t end = clock();
    double time = double(end - start) / CLOCKS_PER_SEC;
    cout << "Time: " << time << " seconds" << endl;

    return 0;
}
