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

#pragma once
#include <map>
#include <vector>
#include <set>
#include <string>
#include "util/io/serializer.h"
#include "util/io/deserializer.h"
#include "basic/value.h"
#include "util/util.h"

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
	static int predefined(const char* s) {
		auto it = rank_map.find(s);
		return it == rank_map.end() ? -1 : (int)it->second;
	}
	static const char* names[count];
	static const std::map<char, char> legacy_rank_codes;
private:
	char r;
	static const std::map<std::string, Rank> rank_map;
	static std::map<std::string, Rank> init_map();
};

struct TaxonomyNodes
{

	TaxonomyNodes(const std::string& file_name);
	TaxonomyNodes(Deserializer &in, uint32_t db_build);
	void save(Serializer &out);
	unsigned get_parent(TaxId taxid) const
	{
		if (taxid < 0 || safe_cast<size_t>(taxid) >= parent_.size())
			return 0;
		return parent_[taxid];
	}
	int rank(TaxId tax_id) const {
		if (tax_id < 0 || safe_cast<size_t>(tax_id) >= rank_.size())
			return -1;
		return rank_[tax_id];
	}
	bool contained(TaxId query, const std::set<TaxId> &filter);
	bool contained(const std::vector<TaxId>& query, const std::set<TaxId> &filter);
	TaxId max() const {
		return safe_cast<TaxId>(parent_.size() - 1);
	}

private:

	uint64_t node_count_;
	std::vector<TaxId> parent_;
	std::vector<Rank> rank_;

};