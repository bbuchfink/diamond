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

#include "tsv.h"
#include "../io/text_input_file.h"

using std::vector;
using std::string;
using std::numeric_limits;
using std::runtime_error;

namespace Util { namespace Tsv {

int64_t count_lines(const std::string& file_name) {
	TextInputFile f(file_name);
	int64_t n = 0;
	while (f.getline(), !f.line.empty() || !f.eof())
		++n;
	f.close();
	return n;
}

}}