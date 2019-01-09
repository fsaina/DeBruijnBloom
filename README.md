Faculty of Electrical Engineering and Computing
Bioinformatics course projects (https://www.fer.unizg.hr/en/course/bio)

Replication of paper:
Space-efficient and exact de Bruijn graph representation based on a Bloom filter ( https://almob.biomedcentral.com/articles/10.1186/1748-7188-8-22 )

### Setup project:
1. Clone this repo
2. Download or compile Jellyfish executable version 2.x.x ( https://github.com/gmarcais/Jellyfish/releases ) into project root directory as 'jellyfish'
3. Compile the project with `g++ -std=c++14 main.cpp BloomFilter.cpp ExactDeBruijnGraph.cpp KmerUtil.cpp MurmurHash3.cpp Tests.cpp -o ./bin/DeBrujinBloom`
4. Run it with `./bin/DeBrujinBloom -k 21 --minAbundance 3 --input data/ecoli/ecoli.fasta`

### Data
The E.Coli genome (data/ecoli.fasta) was taken from http://bacteria.ensembl.org/Escherichia_coli_14a/Info/Index/

To create custom synthetic genome sequences:
1. Make sure `./data/generate.sh` is executable, otherwise run `chmod +x ./data/generate.sh`
2. Run `./data/generate.sh N > ./data/<NAME_OF_FILE>.fasta` (e.g. `./data/generate.sh 100 > ./data/uniform_100.fasta`)

To create the synthetic genome sequences(E.Coli):
1. Compile wgsim tool `gcc -g -O2 -Wall -o ./bin/wgsim ./wgsim/wgsim.c -lz -lm`
2. For creating synthetic sequences of length x, with error y, run
   `./bin/wgsim -1 x -d0 -S11 -e0 -r y ./data/Escherichia_coli.fa ./data/read-x-y.fq /dev/null`

This project is licensed under the terms of the MIT license.
