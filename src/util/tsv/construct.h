#pragma once
#include <stdexcept>
#include "def.h"

namespace Util { namespace Tsv {

template<typename T1, typename... Targs>
struct Construct {
	template<typename T, typename Tok, typename... Uargs>
	T operator()(Tok& tok, Uargs... args) {
		const T1 v = convert_string<T1>(*tok);
		++tok;
		return Construct<Targs...>().template operator()<T, Tok>(tok, args..., v);
	}
};

template<>
struct Construct<void> {
	template<typename T, typename Tok, typename... Uargs>
	T operator()(Tok&, Uargs... args) {
		return T(args...);
	}
};

}}