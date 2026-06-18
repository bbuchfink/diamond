/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <string>
#include "util/io/file.h"
#include "util/string/tokenizer.h"
#include "basic/value.h"
#include "util/string/string.h"

template<typename F>
inline void read_nodes_dmp(const std::string& file_name, F& f)
{
	File file(file_name, "rb", File::Flags::DETECT_COMPRESSION);
	TaxId taxid, parent;
	std::string rank;
	const char* l;
	while (l = file.getline(), !file.eof() || l[0] != '\0') {
		Util::String::Tokenizer<Util::String::StringDelimiter>(l, Util::String::StringDelimiter("\t|\t")) >> taxid >> parent >> rank;
		f(taxid, parent, rank);
	}
	file.close();
}

template<typename F>
inline void read_names_dmp(const std::string& file_name, F& f) {
	File in(file_name, "rb", File::Flags::DETECT_COMPRESSION);
	std::string name, type;
	TaxId id;
	const char* l;
	while (l = in.getline(), !in.eof() || l[0] != '\0') {
		if (std::string(l).empty())
			continue;
		Util::String::Tokenizer<Util::String::StringDelimiter>(l, Util::String::StringDelimiter("\t|\t")) >> id >> name >> Util::String::Skip() >> type;
		type = rstrip(type, "\t|");
		if (type == "scientific name") {
			f(id, name);
		}
	}
	in.close();
}