/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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
#include "../util/io/output_file.h"
#include "../util/data_structures/compact_array.h"

template<typename Type, typename Cmp>
struct ExternalSorter;

struct TaxonList : public CompactArray<vector<unsigned>>
{
	TaxonList(Deserializer &in, size_t size, size_t data_size);
	static void build(OutputFile &db, ExternalSorter<std::pair<std::string, uint32_t>, std::less<std::pair<std::string, uint32_t>>>& accessions, size_t seqs);
};