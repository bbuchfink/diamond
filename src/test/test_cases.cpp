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
{ "blastp (default)", "blastp -p1" },
{ "blastp (multithreaded)", "blastp -p4" },
{ "blastp (blocked)", "blastp -c1 -b0.00002 -p4" },
{ "blastp (more-sensitive)", "blastp --more-sensitive -c1 -p4" },
{ "blastp (target seqs)", "blastp -k3 -c1 -p4" },
{ "blastp (top)", "blastp --top 10 -p4"},
{ "blastp (evalue)", "blastp -e10000 --more-sensitive -c1 -p4" },
{ "blastp (blosum50)", "blastp --matrix blosum50 -p4"},
{ "blastp (pairwise format)", "blastp -c1 -f0 -p4" },
{ "blastp (XML format)", "blastp -c1 -f xml -p4" },
{ "blastp (PAF format)", "blastp -c1 -f paf -p1" }
};

const uint64_t ref_hashes[] = {
0xd03ed9ee31fc5d64,
0xd03ed9ee31fc5d64,
0x361c617388d5cb55,
0x517e7da1d95ff5e4,
0x938b7662efce5ee4,
0x8ac0057c68f8c239,
0x284d7885b1e7f5dc,
0x12885fdee89cd143,
0xbfa21ea5712d871c,
0x2b0922c2d8a6de72,
0x6fa5019017db2bbd,
};

}