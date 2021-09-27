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
0x23b9a0c9b142f36f,
0x23b9a0c9b142f36f,
0x23b9a0c9b142f36f,
0x61d96d8a6a0b431f,
0xcad92644fa4879de,
0x72bae28250e59c6c,
0xe2cfabe9cbda9f18,
0x61d96d8a6a0b431f,
0x684b53c9ae4e539d,
0xd5b8086a38319bde,
0x3e28c839c154b53d,
0x5b318d0a28d30365,
0x5a1f1e76d0dc554b,
0xad9f05902ed6cfc3,
0x201a627d0d128fd5,
0xe1544936bd7b173a,
0xb0d576cc3aadbd3,
0xa2519e06e3bfa2fd,
0x1edd906017a3ddb2,
0x67b3a14cdd541dc3
};

}