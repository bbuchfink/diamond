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

#include <set>
#include <stdexcept>
#include "taxonomy.h"
#include "data/blastdb/taxdmp.h"
#include "basic/config.h"
#include "util/log_stream.h"
#include "taxonomy_nodes.h"

using std::string;
using std::map;
using std::endl;
using std::set;
using std::vector;
using std::runtime_error;

const char* Rank::names[] = {
	"no rank", "superkingdom", "cellular root", "acellular root", "domain", "realm", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "cohort", "subcohort", "superorder",
	"order", "suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus", "section", "subsection", "series", "species group",
	"species subgroup", "species", "subspecies", "varietas", "forma", "strain", "biotype", "clade", "forma specialis", "genotype", "isolate", "morph", "pathogroup", "serogroup", "serotype", "subvariety"
};

map<string, Rank> Rank::init_map() {
	map<string, Rank> r;
	for (size_t i = 0; i < count; ++i)
		r[names[i]] = Rank(i);
	return r;
}

const map<string, Rank> Rank::rank_map = Rank::init_map();

Rank::Rank(const char *s) {
	if (rank_map.find(s) == rank_map.end())
		throw runtime_error("Invalid taxonomic rank: " + string(s));
	r = rank_map.find(s)->second.r;
}

Taxonomy taxonomy;

size_t Taxonomy::load_names() {
	size_t n;
	auto f = [&](int64_t id, const string& name) {
		name_.resize(id + 1);
		name_[id] = name;
		++n;
		};
	read_names_dmp(config.namesdmp, f);
	return n;
}

void Taxonomy::init()
{
	TaskTimer timer;
	if (!config.namesdmp.empty()) {
		timer.go("Loading taxonomy names");
		size_t n = load_names();
		timer.finish();
		message_stream << "Loaded taxonomy names for " << n << " taxon ids." << endl;
	}
}