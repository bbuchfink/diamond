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
0xa941ea1bcaae9cb3,
0xa941ea1bcaae9cb3,
0x7992486f9bc878e8,
0x6f1103a94ffd1a2b,
0xb5b52a7733b476a2,
0xeaf960fbf7c63ff9,
0xa839eaaf7c454ff2,
0x6f1103a94ffd1a2b,
0x4b7ca89df4038f3d,
0x8b983cbaaff963ff,
0x113e1df71d47ee35,
0xc955cd40c085e64d,
0x92297cfae3e80486,
0x563e4f33df3c673d,
0x6a9d5bf640fc1f4b,
0xfc5ca0d04d30faca,
0xf3dadda955ca2d30,
0x45e4056064e260c6,
0xdffb0103534fe08f,
0x778a9e9e5f7a6d64,
};

}