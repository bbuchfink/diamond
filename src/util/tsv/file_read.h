/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

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