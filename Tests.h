#ifndef DEBRUIJNBLOOM_TESTS_H
#define DEBRUIJNBLOOM_TESTS_H

#include <string>
#include <vector>

using namespace std;

namespace Tests {

    void run_all_tests(string inputPath, int k);

    bool create_test_bloom_implementation(string inputPath, int k);

    bool create_test_de_bruijn_graph(string inputPath, int k);

    bool create_kmer_extensions_test();

    bool custom_test_bloom_implementation();
}


#endif
