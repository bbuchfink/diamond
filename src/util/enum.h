#pragma once
#include <map>
#include <string>
#include <stdexcept>

template<typename _t>
using EMap = std::map<_t, std::string>;

template<typename _t> struct EnumTraits {};

template<typename _t> std::string to_string(_t v) {
	auto it = EnumTraits<_t>::to_string.find(v);
	if (it == EnumTraits<_t>::to_string.end())
		throw std::runtime_error("Invalid conversion from enum to string.");
	return it->second;
}
