Faculty of Electrical Engineering and Computing
Bioinformatics projects

Replication of paper:
Space-efficient and exact de Bruijn graph representation based on a Bloom filter ( https://almob.biomedcentral.com/articles/10.1186/1748-7188-8-22 )

### Setup project:
1. Clone this repo
2. Download or compile Jellyfish executable ( https://github.com/gmarcais/Jellyfish/releases ) into project root directory as 'jellyfish'
3. Compile the project with `g++ -std=c++14 main.cpp -o DeBrujinBloom`

### Data
The E.Coli genome (data/ecoli.fasta) was taken from http://bacteria.ensembl.org/Escherichia_coli_14a/Info/Index/
// TODO create synthetic genome sequences