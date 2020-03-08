#include <string.h>
#include <assert.h>
#include "value.h"
#include "../util/util.h"

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