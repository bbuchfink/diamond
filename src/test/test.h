#ifndef TEST_H_
#define TEST_H_

#include <vector>
#include <random>
#include "../basic/sequence.h"

namespace Test {

std::vector<char> generate_random_seq(size_t length, std::minstd_rand0 &rand_engine);
std::vector<char> simulate_homolog(const sequence &seq, double id, std::minstd_rand0 &random_engine);

extern const char* seqs[377][2];

}

#endif