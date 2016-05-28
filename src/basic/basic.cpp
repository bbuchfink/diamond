#include "value.h"
#include "reduction.h"
#include "shape_config.h"
#include "translate.h"
#include "statistics.h"

const char* Const::version_string = "0.8.3";
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
