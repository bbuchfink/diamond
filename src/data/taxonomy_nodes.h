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

#ifndef TAXONOMY_NODES_H_
#define TAXONOMY_NODES_H_

#include <vector>
#include <set>
#include "../util/io/serializer.h"
#include "../util/io/deserializer.h"

using std::vector;
using std::set;

struct TaxonomyNodes
{

	TaxonomyNodes(Deserializer &in);
	static void build(Serializer &out);
	unsigned get_parent(unsigned taxid) const
	{
		if (taxid >= parent_.size())
			throw std::runtime_error(string("No taxonomy node found for taxon id ") + to_string(taxid));
		return parent_[taxid];
	}
	unsigned get_lca(unsigned t1, unsigned t2) const;
	bool contained(unsigned query, const set<unsigned> &filter);
	bool contained(const vector<unsigned> query, const set<unsigned> &filter);

private:

	void set_cached(unsigned taxon_id, bool contained)
	{
		cached_[taxon_id] = true;
		contained_[taxon_id] = contained;
	}

	vector<unsigned> parent_;
	vector<bool> cached_, contained_;

};

#endif