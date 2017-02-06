/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include "value.h"
#include "reduction.h"
#include "shape_config.h"
#include "translate.h"
#include "statistics.h"
#include "sequence.h"

const char* Const::version_string = "0.8.36";
const char* Const::program_name = "diamond";
const char* Const::id_delimiters = " \a\b\f\n\r\t\v";

Value_traits::Value_traits(const char *alphabet, Letter mask_char, const char *ignore) :
	alphabet(alphabet),
	alphabet_size((unsigned)strlen(alphabet)),
	mask_char(mask_char),
	from_char(Char_representation((unsigned)alphabet_size, alphabet, mask_char, ignore))
{}

const char Char_representation::invalid = '\xff';

const Value_traits amino_acid_traits("ARNDCQEGHILKMFPSTWYVBJZX*", 23, "UO-");
const Value_traits nucleotide_traits("ACGTN", 4, "MRWSYKVHDBX");
Value_traits value_traits(amino_acid_traits);
Value_traits input_value_traits(amino_acid_traits);

Align_mode::Align_mode(unsigned mode) :
	mode(mode)
{
	sequence_type = amino_acid;
	switch (mode) {
	case blastx:
		input_sequence_type = nucleotide;
		query_contexts = 6;
		query_translated = true;
		query_len_factor = 3;
		break;
	default:
		input_sequence_type = amino_acid;
		query_contexts = 1;
		query_translated = false;
		query_len_factor = 1;
	}
}

unsigned Align_mode::from_command(unsigned command)
{
	switch (command) {
	case Config::blastx:
		return blastx;
	default:
		return blastp;
	}
}

Align_mode align_mode (Align_mode::blastp);

Reduction Reduction::reduction("KREDQN C G H M F Y ILV W P STA");
//Reduction Reduction::reduction("A KR EDNQ C G H ILVM FYW P ST"); // murphy.10
//const Reduction Reduction::reduction("G D N AEFIKLMQRVW Y H C T S P"); // gmbr.10
//const Reduction Reduction::reduction("EKQR IV LY F AM W HT C DNS"); // dssp.10
//const Reduction Reduction::reduction("K R E D Q N C G H M F Y I L V W P S T A");

Statistics statistics;

const char* shape_codes[][Const::max_shapes] = {
	{ "111101011101111", "111011001100101111", "1111001001010001001111", "111100101000010010010111", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },				// 4x12
	{ "1111011111",		// 16x9
	"111001101111",
	"11101100101011",
	"11010010111011",
	"111010100001111",
	"1110100011001011",
	"11100010100101011",
	"11011000001100111",
	"1101010010000010111",
	"11100001000100100111",
	"110110000100010001101",
	"1110000100001000101011",
	"1101010000010001001011",
	"1101001001000010000111",
	"1101000100100000100000111",
	"1110001000100000001010011" },
	{
		"11001011",		// 16x5
		"101010011",
	"100110101",
	"1110000101",
	"110000100011",
	"1010010000011",
	"1100000010011",
	"11010000000101",
	"100100010000101",
	"1010000000000100011",
	"1010000001000001001",
	"1100000000100001001",
	"10100010000000100001",
	"10010001000000000101",
	"110000000100000010001",
	"10010000100000000000011" },

	{
		"11101011", // 16x6
"110100111",
"11001000111",
"1100001001011",
"10101000010011",
"101001000001011",
"1100010000001011",
"11010000010001001",
"100100100000010101",
"101001000100000101",
"1010001000010000101",
"11001000000100000011",
"101000001000000010011",
"1100010000000100000101",
"11000001000000000100011",
"101000010000000000010011"
	},

	{ "1110010111", // 16x7
	"11001101011",
	"1101001000111",
	"11100010010011",
	"110100101000011",
	"1100100010010101",
	"1101010000010011",
	"1100100000101011",
	"11010001000010011",
	"10101000010001011",
	"11000010010000111",
	"11100000001000001011",
	"110000100010000001101",
	"11010000100000000010011",
	"10100010000010000001011",
	"110001000000010001000101"

	},
	{
		"101011",		// 16x4
		"110011",
	"110000101",
	"1001000011",
	"10010000011",
	"110000010001",
	"1100000001001",
	"10001000000101",
	"10100000100001",
	"100100000000011",
	"101000000010001",
	"1010000001000001",
	"1000010000001001",
	"101000000000000011",
	"100010000000000000101",
	"1000100000000000100001"
	}
};

shape_config shapes;
unsigned shape_from, shape_to;

const Letter Translator::reverseLetter[5] = { 3, 2, 1, 0, 4 };

Letter Translator::lookup[5][5][5];
Letter Translator::lookupReverse[5][5][5];

const Letter Translator::STOP(value_traits.from_char('*'));

const char* Translator::codes[] = {
	0,
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 1
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", // 2
	"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 3
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 4
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", // 5
	"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 6
	0,
	0,
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", // 9
	"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 10
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 11
	"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 12
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", // 13
	"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", // 14
	0,
	"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 16
	0,
	0,
	0,
	0,
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", // 21
	"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 22
	"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 23
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", // 24
	"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 25
	"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" // 26
};

void Translator::init(unsigned id)
{
	static const unsigned idx[] = { 2, 1, 3, 0 };
	if (id >= sizeof(codes) / sizeof(codes[0]) || codes[id] == 0)
		throw std::runtime_error("Invalid genetic code id.");
	for (unsigned i = 0; i < 5; ++i)
		for (unsigned j = 0; j < 5; ++j)
			for (unsigned k = 0; k < 5; ++k)
				if (i == 4 || j == 4 || k == 4) {
					lookup[i][j][k] = value_traits.mask_char;
					lookupReverse[i][j][k] = value_traits.mask_char;
				}
				else {
					lookup[i][j][k] = value_traits.from_char(codes[id][(int)idx[i] * 16 + (int)idx[j] * 4 + (int)idx[k]]);
					lookupReverse[i][j][k] = value_traits.from_char(codes[id][idx[(int)reverseLetter[i]] * 16 + idx[(int)reverseLetter[j]] * 4 + idx[(int)reverseLetter[k]]]);
				}
	for (unsigned i = 0; i < 4; ++i)
		for (unsigned j = 0; j < 4; ++j) {
			if (equal(lookup[i][j], 4))
				lookup[i][j][4] = lookup[i][j][0];
			if (equal(lookupReverse[i][j], 4))
				lookupReverse[i][j][4] = lookupReverse[i][j][0];
		}
}

vector<Letter> sequence::from_string(const char* str)
{
	vector<Letter> seq;
	while (*str)
		seq.push_back(value_traits.from_char(*(str++)));
	return seq;
}

void Seed::enum_neighborhood(unsigned pos, int treshold, vector<Seed>& out, int score)
{
	Letter l = data_[pos];
	score -= score_matrix(l, l);
	for (unsigned i = 0; i < 20; ++i) {
		int new_score = score + score_matrix(l, i);
		data_[pos] = i;
		if (new_score >= treshold) {
			if (pos < config.seed_weight - 1)
				enum_neighborhood(pos + 1, treshold, out, new_score);
			else
				out.push_back(*this);
		}
	}
	data_[pos] = l;
}

void Seed::enum_neighborhood(int treshold, vector<Seed>& out)
{
	out.clear();
	enum_neighborhood(0, treshold, out, score(*this));
}