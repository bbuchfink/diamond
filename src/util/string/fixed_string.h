#pragma once
#include <array>
#include "../algo/MurmurHash3.h"

extern std::array<char, 16> fixed_string_seed;

template<size_t L>
struct FixedString {
	FixedString(const std::string& s) {
		if (s.length() >= L)
			throw std::runtime_error("FixedString");
		std::copy(s.begin(), s.end(), chars.begin());
		chars[s.length()] = '\0';
	}
	bool operator==(const FixedString& s) const {
		return strcmp(chars.data(), s.chars.data()) == 0;
	}
	std::array<char, L> chars;
	struct Hash {
		size_t operator()(const FixedString& s) const {
			std::array<char, 16> out;
			MurmurHash3_x64_128(s.chars.data(), (int)strlen(s.chars.data()), fixed_string_seed.data(), out.data());
			return *(size_t*)out.data();
			//return std::hash<string>()(string(s.chars.data()));
		}
	};
};
