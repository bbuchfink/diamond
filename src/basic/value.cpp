#include <string.h>
#include <assert.h>
#include "value.h"
#include "../util/util.h"

const Letter Char_representation::invalid = '\xff';

invalid_sequence_char_exception::invalid_sequence_char_exception(char ch) :
	msg(std::string("Invalid character (") + print_char(ch) + ") in sequence")
{ }

Char_representation::Char_representation(unsigned size, const char* chars, char mask, const char* mask_chars)
{
	memset(data_, invalid, sizeof(data_));
	for (unsigned i = 0; i < size; ++i) {
		assert(chars[i] != (char)invalid);
		data_[(long)chars[i]] = i;
		data_[(long)tolower(chars[i])] = i;
	}
	while (*mask_chars != 0) {
		const char ch = *mask_chars;
		data_[(long)ch] = mask;
		data_[(long)tolower(ch)] = mask;
		++mask_chars;
	}
}

Value_traits::Value_traits(const char* alphabet, Letter mask_char, const char* ignore) :
	alphabet(alphabet),
	alphabet_size((unsigned)strlen(alphabet)),
	mask_char(mask_char),
	from_char(Char_representation((unsigned)alphabet_size, alphabet, mask_char, ignore))
{}

const Value_traits amino_acid_traits(AMINO_ACID_ALPHABET, 23, "UO-");
const Value_traits nucleotide_traits("ACGTN", 4, "MRWSYKVHDBX");
Value_traits value_traits(amino_acid_traits);
Value_traits input_value_traits(amino_acid_traits);