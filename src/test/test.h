/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once
#include <vector>
#include <random>
#include <stdint.h>
#include <string>
#include <utility>
#include "../basic/sequence.h"

namespace Test {

struct TestCase {
	const char *desc, *command_line;
};

std::vector<Letter> generate_random_seq(size_t length, std::minstd_rand0 &rand_engine);
std::vector<Letter> simulate_homolog(const Sequence &seq, double id, std::minstd_rand0 &random_engine);

extern const std::vector<std::pair<std::string, std::string>> seqs;
extern const std::vector<TestCase> test_cases;
extern const std::vector<uint64_t> ref_hashes;

}