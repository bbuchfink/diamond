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

#include "test.h"

namespace Test {

const TestCase test_cases[] = {
{ "blastp (default)", "blastp" },
{ "blastp (blocked)", "blastp -c1 -b0.00002" },
{ "blastp (more-sensitive, target seqs, evalue)", "blastp -k3 -e10000 --more-sensitive -c1" },
{ "blastp (top, blosum50)", "blastp --top 10 --matrix blosum50"},
{ "blastp (pairwise format)", "blastp -c1 -f0" },
{ "blastp (XML format)", "blastp -c1 -f xml" },
{ "blastp (DAA format)", "blastp -c1 -f daa -p1" }
};

const uint64_t ref_hashes[] = {
0x9bbab2d4498bdaa2,
0x677b5dca34af6ff2,
0x18c00a55f83df537,
0xa1e3026066226494,
0x8a869bfae83a7d1d,
0xdb420eeb9c01821d,
0x834b5c506f0a2c5a
};

}