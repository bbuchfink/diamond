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

std::string get_accession(const string &t);
std::vector<std::string> accession_from_title(const char *title);

struct Taxonomy
{
	void init();
	void load_nodes();
	size_t load_names();

	unsigned get_parent(unsigned taxid) const
	{
		if (taxid >= parent_.size())
			throw std::runtime_error(std::string("No taxonomy node found for taxon id ") + std::to_string(taxid));
		return parent_[taxid];
	}

	unsigned get_lca(unsigned t1, unsigned t2) const;
	
	std::vector<unsigned> parent_;
	std::vector<std::string> name_;
	std::vector<Rank> rank_;

	friend struct TaxonomyNodes;

};

extern Taxonomy taxonomy;