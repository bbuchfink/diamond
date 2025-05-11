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
#include "util/io/serializer.h"
#include "util/io/deserializer.h"
#include "basic/value.h"

struct Rank {
	Rank() :
		r(none)
	{}
	Rank(size_t i) :
		r((char)i)
	{}
	Rank(const char *s);
	enum {
		count = 49, strain = 38, biotype = 39, clade = 40, forma_specialis = 41, genotype = 42, isolate = 43, morph = 44, pathogroup = 45, serogroup = 46, serotype = 47, subvariety = 48,
		forma = 37, varietas = 36, subspecies = 35, species = 34, species_subgroup = 33, species_group = 32, series = 31, subsection = 30, section = 29, subgenus = 28, genus = 27, subtribe = 26,
		tribe = 25, subfamily = 24, family = 23, superfamily = 22, parvorder = 21, infraorder = 20, suborder = 19, order = 18, superorder = 17, subcohort = 16, cohort = 15, infraclass = 14,
		subclass = 13, class_rank = 12, superclass = 11, subphylum = 10, phylum = 9, superphylum = 8, subkingdom = 7, kingdom = 6, realm = 5, domain = 4, acellular_root = 3, cellular_root = 2, superkingdom = 1, none = 0
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
	std::vector<TaxId> lineage(TaxId taxid) const;

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
