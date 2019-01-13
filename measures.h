#ifndef DEBRUIJNBLOOM_MEASURES_H
#define DEBRUIJNBLOOM_MEASURES_H

#include <vector>
#include <string>

using namespace std;

namespace measures {
    float n50(vector<int> lengths);
    vector<int> read_seq_counts(string inputPath);
};


#endif //DEBRUIJNBLOOM_MEASURES_H
