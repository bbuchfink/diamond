/****
DIAMOND protein aligner
Copyright (C) 2024 Max Planck Society for the Advancement of Science e.V.
				   Benjamin Buchfink

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
#include "file.h"
#include "construct.h"

namespace Util { namespace Tsv {

template<typename... Targs, typename Out>
void File::read(Out out) {
	using T = typename Out::container_type::value_type;
	std::unique_ptr<Util::String::TokenizerBase> tok(config_.line_tokenizer->clone());
	if ((int)schema_.size() != sizeof...(Targs))
		throw std::runtime_error("Template parameters do not match schema.");
	rewind();
	while (file_->getline(), !file_->line.empty() || !file_->eof()) {
		tok->reset(file_->line.c_str(), file_->line.c_str() + file_->line.length());
		*out++ = Construct<Targs..., void>().template operator()<T>(tok.get());
	}
}

}}