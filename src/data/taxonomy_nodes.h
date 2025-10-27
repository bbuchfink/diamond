/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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
		if (safe_cast<size_t>(taxid) >= parent_.size())
			throw std::runtime_error(std::string("No taxonomy node found for taxon id ") + std::to_string(taxid));
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

	std::vector<TaxId> parent_;
	std::vector<Rank> rank_;

};