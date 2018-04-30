/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef TAXONOMY_FILTER_H_
#define TAXONOMY_FILTER_H_

#include <set>
#include <vector>
#include <string>
#include "taxon_list.h"
#include "taxonomy_nodes.h"
#include "../util/util.h"

using std::vector;
using std::string;
using std::set;

struct TaxonomyFilter : public vector<bool>
{
	TaxonomyFilter(const string &filter, const TaxonList &list, TaxonomyNodes &nodes)
	{
		const set<unsigned> taxon_filter_list(parse_csv(filter));
		if (taxon_filter_list.empty())
			throw std::runtime_error("Option --taxonlist used with empty list.");
		if (taxon_filter_list.find(1) != taxon_filter_list.end() || taxon_filter_list.find(0) != taxon_filter_list.end())
			throw std::runtime_error("Option --taxonlist used with invalid argument (0 or 1).");
		for (size_t i = 0; i < list.size(); ++i)
			push_back(nodes.contained(list[i], taxon_filter_list));
	}
};

#endif