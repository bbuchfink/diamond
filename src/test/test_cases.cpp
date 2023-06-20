/****
DIAMOND protein aligner
Copyright (C) 2019-2020 Max Planck Society for the Advancement of Science e.V.

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

#include "test.h"

using std::vector;

namespace Test {

const vector<TestCase> test_cases = {
{ "blastp (default)", "blastp -p1" },
{ "blastp (multithreaded)", "blastp -p4" },
{ "blastp (blocked)", "blastp -c1 -b0.00002 -p4" },
{ "blastp (more-sensitive)", "blastp --more-sensitive -c1 -p4" },
{ "blastp (very-sensitive)", "blastp --very-sensitive -c1 -p4" },
{ "blastp (ultra-sensitive)", "blastp --ultra-sensitive -c1 -p4" },
{ "blastp (max-hsps)", "blastp --more-sensitive -c1 -p4 --max-hsps 0" },
{ "blastp (target-parallel)", "blastp --more-sensitive -c1 -p4 --query-parallel-limit 1" },
{ "blastp (query-indexed)", "blastp --more-sensitive -c1 -p4 --algo 1" },
{ "blastp (comp-based-stats 0)", "blastp --more-sensitive -c1 -p4 --comp-based-stats 0" },
{ "blastp (comp-based-stats 2)", "blastp --more-sensitive -c1 -p4 --comp-based-stats 2" },
{ "blastp (comp-based-stats 3)", "blastp --more-sensitive -c1 -p4 --comp-based-stats 3" },
{ "blastp (comp-based-stats 4)", "blastp --more-sensitive -c1 -p4 --comp-based-stats 4" },
{ "blastp (target seqs)", "blastp -k3 -c1 -p4" },
{ "blastp (top)", "blastp --top 10 -p4"},
{ "blastp (evalue)", "blastp -e10000 --more-sensitive -c1 -p4" },
{ "blastp (blosum50)", "blastp --matrix blosum50 -p4"},
{ "blastp (pairwise format)", "blastp -c1 -f0 -p4" },
{ "blastp (XML format)", "blastp -c1 -f xml -p4" },
{ "blastp (PAF format)", "blastp -c1 -f paf -p1" }
};

const vector<uint64_t> ref_hashes = {
0x36bf16afef49c7ad,
0x36bf16afef49c7ad,
0x36bf16afef49c7ad,
0x7ed13391c638dc2e,
0x61ac7ee1bb73d36d,
0xd62b1c97fb27608f,
0x2dd4b2985c1bebd2,
0x7ed13391c638dc2e,
0x7ed13391c638dc2e,
0x9a20976998759371,
0xa67de9d0530d5968,
0xa67de9d0530d5968,
0x3d593e440ca8eb97,
0x487a213a131d4958,
0x201a627d0d128fd5,
0xe787dcb23cc5b120,
0x5aa4baf48a888be9,
0xa2519e06e3bfa2fd,
0xe65d44d7b1055824,
0x67b3a14cdd541dc3,
};

}