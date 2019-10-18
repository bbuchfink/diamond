/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef METADATA_H_
#define METADATA_H_

#include <vector>
#include <string>
#include "taxon_list.h"
#include "taxonomy_nodes.h"
#include "taxonomy.h"

struct Metadata
{
	Metadata():
		taxon_list(nullptr),
		taxon_nodes(nullptr),
		taxon_filter(nullptr),
		taxonomy_scientific_names(nullptr)
	{}
	void free()
	{
		delete taxon_list;
		delete taxon_nodes;
		delete taxon_filter;
		delete taxonomy_scientific_names;
		taxon_list = nullptr;
		taxon_nodes = nullptr;
		taxon_filter = nullptr;
		taxonomy_scientific_names = nullptr;
	}
	TaxonList *taxon_list;
	TaxonomyNodes *taxon_nodes;
	TaxonomyFilter *taxon_filter;
	std::vector<std::string> *taxonomy_scientific_names;
};

#endif