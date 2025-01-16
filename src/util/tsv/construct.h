/****
DIAMOND protein aligner
Copyright (C) 2019-2024 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#pragma once
#include "util/string/tokenizer_dyn.h"

namespace Util { namespace Tsv {

template<typename T1, typename... Targs>
struct Construct {
	template<typename T, typename... Uargs>
	T operator()(Util::String::TokenizerBase* tok, Uargs... args) {
		const T1 v = Util::String::convert_string<T1>(**tok);
		++(*tok);
		return Construct<Targs...>().template operator()<T>(tok, args..., v);
	}
};

template<>
struct Construct<void> {
	template<typename T, typename... Uargs>
	T operator()(Util::String::TokenizerBase*, Uargs... args) {
		return T(args...);
	}
};

}}