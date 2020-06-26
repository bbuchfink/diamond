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
{ "blastp (target-parallel)", "blastp --more-sensitive -c1 -p4 --query-parallel-limit 1" },
{ "blastp (query-indexed)", "blastp --more-sensitive -c1 -p4 --algo 1" },
{ "blastp (target seqs)", "blastp -k3 -c1 -p4" },
{ "blastp (top)", "blastp --top 10 -p4"},
{ "blastp (evalue)", "blastp -e10000 --more-sensitive -c1 -p4" },
{ "blastp (blosum50)", "blastp --matrix blosum50 -p4"},
{ "blastp (pairwise format)", "blastp -c1 -f0 -p4" },
{ "blastp (XML format)", "blastp -c1 -f xml -p4" },
{ "blastp (PAF format)", "blastp -c1 -f paf -p1" }
};

const vector<uint64_t> ref_hashes = {
0x4f959b7506a0b621,
0x4f959b7506a0b621,
0x2d95a36e0a13bf9b,
0x959b05442e7bfebf,
0x959b05442e7bfebf,
0x3bd7a39c0bd45a3,
0x8ac0057c68f8c239,
0xcc42d6752cbc88c1,
0xa76a583e6fe4c891,
0xc40f9bffc00c4f97,
0x1d42766969a3525d,
0xdf32d46bc18e400a,
};

}