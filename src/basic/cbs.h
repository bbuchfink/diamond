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
    double lambda_ratio;

};