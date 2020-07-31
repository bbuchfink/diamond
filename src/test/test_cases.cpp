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
{ "blastp (comp-based-stats)", "blastp --more-sensitive -c1 -p4 --comp-based-stats 0" },
{ "blastp (target seqs)", "blastp -k3 -c1 -p4" },
{ "blastp (top)", "blastp --top 10 -p4"},
{ "blastp (evalue)", "blastp -e10000 --more-sensitive -c1 -p4" },
{ "blastp (blosum50)", "blastp --matrix blosum50 -p4"},
{ "blastp (pairwise format)", "blastp -c1 -f0 -p4" },
{ "blastp (XML format)", "blastp -c1 -f xml -p4" },
{ "blastp (PAF format)", "blastp -c1 -f paf -p1" }
};

const vector<uint64_t> ref_hashes = {
0xea2c4d83866391cc,
0xea2c4d83866391cc,
0x5534fb26c4951ea8,
0xb08613ba4d95d8e5,
0x322ccc108a8fcd4f,
0xd248fd6f8a681713,
0xa682a0031cf6e942,
0xb08613ba4d95d8e5,
0x3b287684ae0261e,
0xe16c417914681afd,
0x8ac0057c68f8c239,
0xb06d5eaf6438c123,
0x3953a4a82ee229b8,
0x4239facd8427148a,
0xa22dad750780029f,
0xb51769cee7147e17
};

}