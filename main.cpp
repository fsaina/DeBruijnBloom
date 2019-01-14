#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

#include "cmdline.h"
#include "ExactDeBruijnGraph.h"
#include "measures.h"


using namespace std;

/*
 * Read a file containing the kmers. Method skips lines that contain mer
 * occurrence count.
 */
vector<string> read_mers(string inputPath){
    ifstream in(inputPath);
    vector<string> mers;

    string line;
    while(getline(in, line)) {
        if (line.find(">") != 0) {
           mers.push_back(line);
        }
    }

    in.close();
    return mers;
}

/*
 * Entry point of the command line application.
 * For parameters required to run this CLI application please reffer to the
 * proved README.md or run the tool with a --help flag.
 */
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
    p.add<bool>("n50", 'n', "Output the N50 measure at the end of execution", false, false);

    p.parse_check(argc, argv);

    int k = p.get<int>("kmers");
    int minAbundance = p.get<int>("minAbundance");
    int maxBreadth = p.get<int>("maxBreadth");
    int maxDepth = p.get<int>("maxDepth");
    string inputPath = p.get<string>("input");
    string outputPath = p.get<string>("output");
    string jellyfishBinPath = p.get<string>("jellyfish");
    string tmpDir = p.get<string>("tmp");
    bool outputN50 = p.get<bool>("n50");

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
              + " -s 100M -t 4 -L " + to_string(minAbundance)
              + " " + inputPath
              + " -o " + tmpDir + '/' + defaultJellyfishOutput;
    system(command.c_str());

    // Convert binary jellyfish output to human readable format
    command = jellyfishBinPath + " dump " + tmpDir + '/' + defaultJellyfishOutput + " > " + tmpDir + '/' + jellyfishTmpFileName;
    system(command.c_str());

    // read the file and load only the k-mers (not their counts)
    inputPath = tmpDir + '/' + jellyfishTmpFileName;
    vector<string> kmers = read_mers(inputPath);

    cout<< "Number of kmers: " << kmers.size() << endl;

    command = "rm -rf " + outputPath;
    system(command.c_str());

    command = "mkdir " + outputPath;
    system(command.c_str());

    string fullOutputPath = outputPath + "/" + defaultProgramOutput;

    ExactDeBruijnGraph graph = ExactDeBruijnGraph(inputPath, kmers.size(), k);
    graph.traverse(inputPath, fullOutputPath, maxBreadth, maxDepth);

    clock_t end = clock();
    double time = double(end - start) / CLOCKS_PER_SEC;
    cout << "Time: " << time << " seconds" << endl;
    cout << "DeBruijn graph size(bytes): " << graph.graphSizeInBytes() << endl;

    if (outputN50 == true) {
        cout << endl;
        vector<int> lengths = measures::read_seq_counts(fullOutputPath);
        float n50_measure = measures::n50(lengths);
        cout << "N50: " << n50_measure << endl;
    }

    return 0;
}
