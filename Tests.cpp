#include <iostream>
#include <cstdio>
#include "BloomFilter.h"
#include "ExactDeBruijnGraph.h"
#include "KmerUtil.h"
#include "Tests.h"


bool Tests::custom_test_bloom_implementation() {
    BloomFilter bf = BloomFilter(42, 6);

    bf.add("abc");
    bf.add("bcd");
    bf.add("cde");
    bf.add("def");

    bool contains;
    // NOTE in bloom filter data structure false positives are possible, but false negatives should NEVER happen

    contains = bf.contains("abc");
    cout << "Should be true : " << boolalpha  << contains << endl;
    contains = bf.contains("bcd");
    cout << "Should be true : " << boolalpha  << contains << endl;
    contains = bf.contains("cde");
    cout << "Should be true : " << boolalpha  << contains << endl;
    contains = bf.contains("def");
    cout << "Should be true : " << boolalpha  << contains << endl;

    contains = bf.contains("cba");
    cout << "Should be false : " << boolalpha << contains << endl;
    contains = bf.contains("dcb");
    cout << "Should be false : " << boolalpha << contains << endl;
    contains = bf.contains("edc");
    cout << "Should be false : " << boolalpha << contains << endl;
    contains = bf.contains("fed");
    cout << "Should be false : " << boolalpha << contains << endl;
    return true;
}

bool Tests::create_test_bloom_implementation(vector<string> &kmers, int k) {
    BloomFilter bf = BloomFilter(kmers, k);

    for (string s : kmers) {
        bf.add(s);
    }

    // NOTE test that there are no false negatives
    for (string s : kmers) {
        bool contains = bf.contains(s);
        if (!contains) {
            printf("> ERROR\n");
        }
    }

    bool contains;
    // NOTE in bloom filter data structure false positives are possible, but false negatives should NEVER happen

    contains = bf.contains("AAAAAAAAAAAAAAAAAAAAA");
    cout << "Should be true : " << boolalpha  << contains << endl;

    contains = bf.contains("ATGAACGGCAGCGGGCCAAAA");
    cout << "Should be false : " << boolalpha << contains << endl;
    return true;
}

bool Tests::create_test_de_bruijn_graph(vector<string> kmers, int k) {
    ExactDeBruijnGraph graph = ExactDeBruijnGraph(kmers, k);
    return true;
}

bool Tests::create_kmer_extensions_test() {
    for (string s : KmerUtil::generateExtensions("ACT")) {
        cout << s << endl;
    }

    return true;
}

void Tests::run_all_tests(vector<string> &kmers, int k) {
    custom_test_bloom_implementation();

    create_test_bloom_implementation(kmers, k);

    create_test_de_bruijn_graph(kmers, k);

    create_kmer_extensions_test();
}
