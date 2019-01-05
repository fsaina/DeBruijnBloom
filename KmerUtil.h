#ifndef DEBRUIJNBLOOM_KMERUTIL_H
#define DEBRUIJNBLOOM_KMERUTIL_H

#include <string>

using namespace std;

namespace KmerUtil {

    const string reverseComplement(string kmer);

    char complement(char c);
}


#endif //DEBRUIJNBLOOM_KMERUTIL_H
