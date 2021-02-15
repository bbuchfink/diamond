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
0x602762c977aa8682,
0x602762c977aa8682,
0x38498d4f4d3eb7c9,
0x44d8f0f470123331,
0x544de0f64c7b7b0,
0x9164580efdd15995,
0x9a54b156f8f2146a,
0x44d8f0f470123331,
0xc5e10fd8002fb6eb,
0x794cb4944f11ffdc,
0xdfe1489c8ea1b4b6,
0xd0350017fe8f8fda,
0xc56f46f150f65fc1,
0xf1274743d0f712bc,
0x5298dd163b9666b3,
0xd20b29b1abecd9c4,
0xae7bc1145b22152f,
0xc43258834622128e,
0x5d81f357aacf347c,
0x58c74e056adf9a71,
};

}