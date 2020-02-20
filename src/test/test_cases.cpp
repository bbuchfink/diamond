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
{ "blastp (default)", "blastp -p4" },
{ "blastp (blocked)", "blastp -c1 -b0.00002 -p4" },
{ "blastp (more-sensitive)", "blastp --more-sensitive -c1 -p4" },
{ "blastp (target seqs)", "blastp -k3 -c1 -p4" },
{ "blastp (top)", "blastp --top 10 -p4"},
{ "blastp (evalue)", "blastp -e10000 --more-sensitive -c1 -p4" },
{ "blastp (blosum50)", "blastp --matrix blosum50 -p4"},
{ "blastp (pairwise format)", "blastp -c1 -f0 -p4" },
{ "blastp (XML format)", "blastp -c1 -f xml -p4" },
{ "blastp (DAA format)", "blastp -c1 -f daa -p1" },
{ "blastp (PAF format)", "blastp -c1 -f paf -p1" }
};

const uint64_t ref_hashes[] = {
0x9bbab2d4498bdaa2,
0x677b5dca34af6ff2,
0xc6bd7cda65df9f6c,
0x36fcb20b94621955,
0x1a606470f4f55e40,
0xbe9037e2769feb1,
0x9bb6322d7466227e,
0x8a869bfae83a7d1d,
0xdb420eeb9c01821d,
0x834b5c506f0a2c5a,
0x9660705d7f35297c
};

}