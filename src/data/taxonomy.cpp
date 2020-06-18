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
#include "../util/merge_sort.h"
#include "../util/log_stream.h"
#include "reference.h"
#include "../util/string/string.h"
#include "../util/string/tokenizer.h"

using namespace std;

const char* Rank::names[] = {
	"no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "cohort", "subcohort", "superorder",
	"order", "suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus", "section", "subsection", "series", "species group",
	"species subgroup", "species", "subspecies", "varietas", "forma", "strain"
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
		throw std::runtime_error("Invalid taxonomic rank.");
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
	return t;
}

void Taxonomy::load()
{
	char acc[max_accesion_len + 2];
	unsigned taxid;
	TextInputFile f(config.prot_accession2taxid);
	f.getline();
	
	while (!f.eof() && (f.getline(), !f.line.empty())) {
		if (sscanf(f.line.c_str(), "%*s%15s%u%*u", acc, &taxid) != 2) {
			//std::cout << f.line << endl;
			throw std::runtime_error("Invalid taxonomy mapping file format.");
		}
		if (strlen(acc) > max_accesion_len) {
			//std::cout << f.line << endl;
			throw std::runtime_error("Accession exceeds supported length.");
		}
		accession2taxid_.push_back(std::make_pair(Accession(acc), taxid));
		/*if (f.line_count % 10000 == 0)
			std::cout << f.line_count << endl;*/
	}
	f.close();
	merge_sort(accession2taxid_.begin(), accession2taxid_.end(), config.threads_);
}

void Taxonomy::load_nodes()
{
	TextInputFile f(config.nodesdmp);
	unsigned taxid, parent;
	string rank;
	while (!f.eof() && (f.getline(), !f.line.empty())) {
		Util::String::Tokenizer(f.line, "\t|\t") >> taxid >> parent >> rank;
		parent_.resize(taxid + 1);
		parent_[taxid] = parent;
		rank_.resize(taxid + 1);
		rank_[taxid] = Rank(rank.c_str());
	}
	f.close();
}

size_t Taxonomy::load_names() {
	TextInputFile in(config.namesdmp);
	string name, type;
	long id;
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
	if (!config.prot_accession2taxid.empty()) {
		timer.go("Loading taxonomy");
		load();
		timer.finish();
		message_stream << "Accession mappings = " << accession2taxid_.size() << endl;
	}
	if (!config.nodesdmp.empty()) {
		timer.go("Loading taxonomy nodes");
		load_nodes();
		timer.finish();
	}
	if (!config.namesdmp.empty()) {
		timer.go("Loading taxonomy names");
		size_t n = load_names();
		timer.finish();
		message_stream << "Loaded taxonomy names for " << n << " taxon ids." << endl;
	}
}

vector<string> Taxonomy::Accession::from_title(const char *title)
{
	vector<string> t(seq_titles(title));
	for (vector<string>::iterator i = t.begin(); i < t.end(); ++i)
		*i = get_accession(blast_id(*i));
	return t;
}

void Taxonomy::get_taxids(const char *id, set<unsigned> &taxons) const
{
	const vector<string> t(seq_titles(id));
	for (vector<string>::const_iterator i = t.begin(); i < t.end(); ++i) {
		const unsigned id = get(Taxonomy::Accession(*i));
		if(id != 0)
			taxons.insert(id);
	}
}

unsigned Taxonomy::get_lca(unsigned t1, unsigned t2) const
{
	static const int max = 64;
	if (t1 == t2 || t2 == 0)
		return t1;
	if (t1 == 0)
		return t2;
	unsigned p = t2;
	set<unsigned> l;
	int n = 0;
	do {
		p = get_parent(p);
		if (p == 0)
			return t1;
		l.insert(p);
		if (++n > max)
			throw std::runtime_error("Path in taxonomy too long (1).");
	} while (p != t1 && p != 1);
	if (p == t1)
		return p;
	p = t1;
	n = 0;
	while (l.find(p) == l.end()) {
		p = get_parent(p);
		if (p == 0)
			return t2;
		if (++n > max)
			throw std::runtime_error("Path in taxonomy too long (2).");
	}
	return p;
}