#pragma once

#include <array>
#include <vector>
#include "sequence.h"

std::array<double, 20> composition(const sequence& s);

struct TargetMatrix {

    TargetMatrix() :
        lambda_ratio(1.0)
    {}

    TargetMatrix(const double* query_comp, const sequence& target);

    std::vector<int8_t> scores;
    std::vector<int32_t> scores32;
    double lambda_ratio;

};

extern const int ALPH_TO_NCBI[20];
extern const double BLOSUM62_FREQRATIOS[28][28];
constexpr double BLOSUM62_UNGAPPED_LAMBDA = 0.3176;