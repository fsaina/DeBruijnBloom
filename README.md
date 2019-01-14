Faculty of Electrical Engineering and Computing
Bioinformatics course projects (https://www.fer.unizg.hr/en/course/bio)

Replication of paper:
Space-efficient and exact de Bruijn graph representation based on a Bloom filter ( https://almob.biomedcentral.com/articles/10.1186/1748-7188-8-22 )

### Setup project
Run `sh ./INSTALL.sh` to setup the project. The script will:
1. Create ./bin directory
2. Download the Jellyfish executable (version 2.2.10 for Linux) into ./bin directory
3. Compile the project with `g++ -std=c++14 main.cpp BloomFilter.cpp ExactDeBruijnGraph.cpp KmerUtil.cpp MurmurHash3.cpp measures.cpp -o ./bin/DeBrujinBloom`

### Run project
Execution example: `./bin/DeBrujinBloom -k 21 --minAbundance 3 --maxBreadth 20 --maxDepth 500 --input data/ecoli/ecoli.fasta`
For more information run `./bin/DeBrujinBloom --help`

### Run test example
To run test example:
1. Run `g++ -std=c++14 Demo.cpp -o bin/Demo` to compile
2. Run `./bin/Demo` to execute

### Data
The E.Coli genome (data/ecoli.fasta) was taken from http://bacteria.ensembl.org/Escherichia_coli_14a/Info/Index/

To create custom synthetic genome sequences:
1. Make sure `./data/generate.sh` is executable, otherwise run `chmod +x ./data/generate.sh`
2. Run `./data/generate.sh N > ./data/<NAME_OF_FILE>.fasta` (e.g. `./data/generate.sh 100 > ./data/uniform_100.fasta`)

To create the synthetic genome sequences(E.Coli):
1. Compile wgsim tool `gcc -g -O2 -Wall -o ./bin/wgsim ./wgsim/wgsim.c -lz -lm`
2. For creating synthetic sequences of length x, with error y, run
   `./bin/wgsim -1 x -d0 -S11 -e0 -r y ./data/Escherichia_coli.fa ./data/read-x-y.fq /dev/null`

### Testing
For testing purposes we used 'blastn' and 'valgrind' tools, and made our implementation of N50 measure(see in next section).

To use blastn(used to check matching of generated sequences and original data):
1. Download and compile blastn and put it in bin/ directory
2. Generate some outputs in fasta format
3. Run `./bin/blastn -query path-to-generated-data.fasta -subject path-to-original-sequence.fasta -out results.txt`

To use Valgrind tool(used to get a memory usage report):
1. Download and install Valgrind
2. Run `valgrind ./bin/DeBrujinBloom program_arguments`

We downloaded the original implementation of this paper at http://minia.genouest.org/, so we can test their performances:
1. Download and compile minia, move it to /bin directory
2. Run `./bin/minia path-to-fasta-file k min-abundance genome-length output-file-name`

### Measure
To calculate the N50 measure from the execution output file, either:
1. Add a "-n 1" flag when calling `./bin/DeBrujinBloom`
2. Or run the N50.py python script with `grep "^>" ./output/output.fasta | cut -d "_" -f5 | ./N50.py` (Don't forget to make N50.py executable with `chmod +x`)

This project is licensed under the terms of the MIT license.
