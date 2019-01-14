#include <fstream>
#include <iostream>
#include <string>

using namespace std;

/*
 * Read a sequence from file. Method skips lines that contain mer
 * occurrence count.
 *
 * Author(s): Marin Vernier
 */
string read_sequence(string inputPath){
    ifstream in(inputPath);
    string sequence;

    string line;
    while(getline(in, line)) {
        if (line.find(">") != 0) {
            sequence = line;
        }
    }

    in.close();
    return sequence;
}

/*
 * Entry point of the demo application.
 *
 * Author(s): Marin Vernier
 */
int main(int argc, char *argv[]) {
    string command = "rm -rf temp";
    system(command.c_str());
    command = "mkdir temp";
    system(command.c_str());

    command = "./bin/DeBrujinBloom -k 21 -a 1 --maxBreadth 20 --maxDepth 500 "
              "-i data/synthetic/synthetic_ecoli-100.fasta -o temp -n 1";
    system(command.c_str());

    string originalSequence = read_sequence("data/synthetic/synthetic_ecoli-100.fasta");
    string generatedSequence = read_sequence("temp/output.fasta");

    cout << "Original: " << originalSequence << endl;
    cout << "Traverse: " << generatedSequence << endl;

    cout << "Are equal? " << boolalpha << (originalSequence.compare(generatedSequence) == 0) << endl;

    return 0;
}
