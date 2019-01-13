wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-macosx -O ./bin/jellyfish
chmod +x ./bin/jellyfish
g++ -std=c++14 main.cpp BloomFilter.cpp ExactDeBruijnGraph.cpp KmerUtil.cpp MurmurHash3.cpp measures.cpp -o ./bin/DeBrujinBloom
