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
#include <map>
#include <vector>
#include <set>
#include <string>
#include "../util/io/serializer.h"
#include "../util/io/deserializer.h"
#include "../basic/value.h"

struct Rank {
	Rank() :
		r(none)
	{}
	Rank(size_t i) :
		r((char)i)
	{}
	Rank(const char *s);
	enum {
		count = 45, strain = 34, biotype = 35, clade = 36, forma_specialis = 37, genotype = 38, isolate = 39, morph = 40, pathogroup = 41, serogroup = 42, serotype = 43, subvariety = 44,
		forma = 33, varietas = 32, subspecies = 31, species = 30, species_subgroup = 29, species_group = 28, series = 27, subsection = 26, section = 25, subgenus = 24, genus = 23, subtribe = 22,
		tribe = 21, subfamily = 20, family = 19, superfamily = 18, parvorder = 17, infraorder = 16, suborder = 15, order = 14, superorder = 13, subcohort = 12, cohort = 11, infraclass = 10,
		subclass = 9, class_rank = 8, superclass = 7, subphylum = 6, phylum = 5, superphylum = 4, subkingdom = 3, kingdom = 2, superkingdom = 1, none = 0
	};
	operator int() const {
		return (int)r;
	}
	friend std::ostream& operator<<(std::ostream &s, Rank &r) {
		s << std::string(names[(int)r.r]);
		return s;
	}
	static const char* names[count];
private:
	char r;
	static const std::map<std::string, Rank> rank_map;
	static std::map<std::string, Rank> init_map();
};

struct TaxonomyNodes
{

	TaxonomyNodes(const std::string& file_name, const bool init_cache = false);
	TaxonomyNodes(Deserializer &in, uint32_t db_build);
	void save(Serializer &out);
	unsigned get_parent(unsigned taxid) const
	{
		if (taxid >= parent_.size())
			throw std::runtime_error(std::string("No taxonomy node found for taxon id ") + std::to_string(taxid));
		return parent_[taxid];
	}
	unsigned rank_taxid(unsigned taxid, Rank rank) const;
	std::set<TaxId> rank_taxid(const std::vector<TaxId> &taxid, Rank rank) const;
	unsigned get_lca(unsigned t1, unsigned t2) const;
	bool contained(TaxId query, const std::set<TaxId> &filter);
	bool contained(const std::vector<TaxId>& query, const std::set<TaxId> &filter);

private:

	void set_cached(unsigned taxon_id, bool contained)
	{
		cached_[taxon_id] = true;
		contained_[taxon_id] = contained;
	}
	void init_cache();

	std::vector<TaxId> parent_;
	std::vector<Rank> rank_;
	std::vector<bool> cached_, contained_;

};
