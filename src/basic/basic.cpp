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

const char* Const::version_string = "0.8.21";
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

const Reduction Reduction::reduction("KREDQN C G H M F Y ILV W P STA");

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
	//{ "11111", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, }
};

shape_config shapes;
unsigned shape_from, shape_to;

const Letter Translator::reverseLetter[5] = { 3, 2, 1, 0, 4 };

const Letter Translator::lookup[5][5][5] = {
	{ { 11,2,11,2,23 },
	{ 16,16,16,16,16 },
	{ 1,15,1,15,23 },
	{ 9,9,12,9,23 },
	{ 23,23,23,23,23 },
	},
	{ { 5,8,5,8,23 },
	{ 14,14,14,14,14 },
	{ 1,1,1,1,1 },
	{ 10,10,10,10,10 },
	{ 23,23,23,23,23 },
	},
	{ { 6,3,6,3,23 },
	{ 0,0,0,0,0 },
	{ 7,7,7,7,7 },
	{ 19,19,19,19,19 },
	{ 23,23,23,23,23 },
	},
	{ { 24,18,24,18,23 },
	{ 15,15,15,15,15 },
	{ 24,4,17,4,23 },
	{ 10,13,10,13,23 },
	{ 23,23,23,23,23 },
	},
	{ { 23,23,23,23,23 },
	{ 23,23,23,23,23 },
	{ 23,23,23,23,23 },
	{ 23,23,23,23,23 },
	{ 23,23,23,23,23 },
	} };

const Letter Translator::lookupReverse[5][5][5] = {
	{ { 13,10,13,10,23 },
	{ 4,17,4,24,23 },
	{ 15,15,15,15,15 },
	{ 18,24,18,24,23 },
	{ 23,23,23,23,23 },
	},
	{ { 19,19,19,19,19 },
	{ 7,7,7,7,7 },
	{ 0,0,0,0,0 },
	{ 3,6,3,6,23 },
	{ 23,23,23,23,23 },
	},
	{ { 10,10,10,10,10 },
	{ 1,1,1,1,1 },
	{ 14,14,14,14,14 },
	{ 8,5,8,5,23 },
	{ 23,23,23,23,23 },
	},
	{ { 9,12,9,9,23 },
	{ 15,1,15,1,23 },
	{ 16,16,16,16,16 },
	{ 2,11,2,11,23 },
	{ 23,23,23,23,23 },
	},
	{ { 23,23,23,23,23 },
	{ 23,23,23,23,23 },
	{ 23,23,23,23,23 },
	{ 23,23,23,23,23 },
	{ 23,23,23,23,23 },
	} };

const Letter Translator::STOP(value_traits.from_char('*'));

vector<Letter> sequence::from_string(const char* str)
{
	vector<Letter> seq;
	while (*str)
		seq.push_back(value_traits.from_char(*(str++)));
	return seq;
}