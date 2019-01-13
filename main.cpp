#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <map>

#include "cmdline.h"
#include "ExactDeBruijnGraph.h"


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
 * Read a fasta file and retrieve the sequence counts (lengths) assigned to
 * each sequence
 */
vector<int> read_seq_counts(string inputPath){
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

/*
 * Returns the N50 value of the passed list of numbers.
 *
 * Based on the Broad Institute definition:
 * https://www.broad.harvard.edu/crd/wiki/index.php/N50
 */
float measure_n50(vector<int> lengths) {
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
    graph.simple_traverse(inputPath, fullOutputPath, maxBreadth, maxDepth);

    clock_t end = clock();
    double time = double(end - start) / CLOCKS_PER_SEC;
    cout << "Time: " << time << " seconds" << endl;

    cout << endl;

    vector<int> lengths = read_seq_counts(fullOutputPath);
    float n50_measure = measure_n50(lengths);
    cout << "N50: " << n50_measure << endl;

    return 0;
}
