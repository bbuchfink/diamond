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

#include <stdio.h>
#include <set>
#include <stdexcept>
#include "taxonomy.h"
#include "../util/io/text_input_file.h"
#include "../basic/config.h"
#include "../util/log_stream.h"
#include "reference.h"
#include "../util/string/string.h"
#include "../util/string/tokenizer.h"
#include "../util/util.h"
#include "../util/sequence/sequence.h"

using std::string;
using std::map;
using std::endl;
using std::set;
using std::vector;

const char* Rank::names[] = {
	"no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "cohort", "subcohort", "superorder",
	"order", "suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus", "section", "subsection", "series", "species group",
	"species subgroup", "species", "subspecies", "varietas", "forma", "strain", "biotype", "clade", "forma specialis", "genotype", "isolate", "morph", "pathogroup", "serogroup", "serotype", "subvariety"
};

map<std::string, Rank> Rank::init_map() {
	map<string, Rank> r;
	for (size_t i = 0; i < count; ++i)
		r[names[i]] = Rank(i);
	return r;
}

const map<std::string, Rank> Rank::rank_map = Rank::init_map();

Rank::Rank(const char *s) {
	if (rank_map.find(s) == rank_map.end())
		throw std::runtime_error("Invalid taxonomic rank: " + string(s));
	r = rank_map.find(s)->second.r;
}

Taxonomy taxonomy;

string get_accession(const string &title)
{
	size_t i;
	string t(title);
	if (t.compare(0, 6, "UniRef") == 0)
		t.erase(0, t.find('_', 0) + 1);
	else if ((i = t.find_first_of('|', 0)) != string::npos) {
		if (t.compare(0, 3, "gi|") == 0) {
			t.erase(0, t.find_first_of('|', i + 1) + 1);
			i = t.find_first_of('|', 0);
		}
		t.erase(0, i + 1);
		i = t.find_first_of('|', 0);
		if (i != string::npos)
			t.erase(i);
	}
	i = t.find_last_of('.');
	if (i != string::npos)
		t.erase(i);
	return t;
}

size_t Taxonomy::load_names() {
	TextInputFile in(config.namesdmp);
	string name, type;
	int64_t id;
	size_t n = 0;
	while (in.getline(), !in.eof()) {
		if (in.line.empty())
			continue;
		Util::String::Tokenizer(in.line, "\t|\t") >> id >> name >> Util::String::Skip() >> type;
		rstrip(type, "\t|");
		if (type == "scientific name") {
			name_.resize(id + 1);
			name_[id] = name;
			++n;
		}
	}
	in.close();
	return n;
}

void Taxonomy::init()
{
	task_timer timer;
	if (!config.namesdmp.empty()) {
		timer.go("Loading taxonomy names");
		size_t n = load_names();
		timer.finish();
		message_stream << "Loaded taxonomy names for " << n << " taxon ids." << endl;
	}
}

vector<string> accession_from_title(const char *title)
{
	vector<string> t(seq_titles(title));
	for (vector<string>::iterator i = t.begin(); i < t.end(); ++i)
		*i = get_accession(Util::Seq::seqid(i->c_str(), false));
	return t;
}
