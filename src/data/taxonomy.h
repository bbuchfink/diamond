/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
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
#include <vector>
#include <algorithm>
#include <string>
#include <set>
#include <map>
#include <ostream>
#include "../basic/const.h"
#include "taxon_list.h"
#include "taxonomy_nodes.h"
#include "../util/data_structures/bit_vector.h"

std::string get_accession(const std::string &t);
std::vector<std::string> accession_from_title(const char *title);

struct Taxonomy
{
	void init();
	size_t load_names();
	
	std::vector<std::string> name_;

	friend struct TaxonomyNodes;

};

extern Taxonomy taxonomy;