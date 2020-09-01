/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include <vector>
#include "../../basic/sequence.h"
#include "../../basic/reduction.h"
#include "../../basic/shape_config.h"

struct SeedTable {

	SeedTable(const sequence& seq, const Reduction& reduction, const shape_config& shapes);

private:

	enum { WORD_SIZE = 6 };

	const size_t size;
	std::vector<uint32_t> lookup_, ptr_;

};