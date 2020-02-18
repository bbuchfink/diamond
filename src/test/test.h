/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef TEST_H_
#define TEST_H_

#include <vector>
#include <random>
#include <stdint.h>
#include "../basic/sequence.h"

namespace Test {

struct TestCase {
	const char *desc, *command_line;
};

std::vector<char> generate_random_seq(size_t length, std::minstd_rand0 &rand_engine);
std::vector<char> simulate_homolog(const sequence &seq, double id, std::minstd_rand0 &random_engine);

enum { CASE_COUNT = 7, PROT_COUNT = 387 };

extern const char* seqs[PROT_COUNT][2];
extern const TestCase test_cases[CASE_COUNT];
extern const uint64_t ref_hashes[CASE_COUNT];

}

#endif