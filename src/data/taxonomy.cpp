/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#include <stdio.h>
#include <set>
#include "taxonomy.h"
#include "../util/compressed_stream.h"
#include "../basic/config.h"
#include "../util/merge_sort.h"
#include "../util/log_stream.h"

using std::cout;
using std::endl;
using std::set;

Taxonomy taxonomy;

string& get_accession(string &t)
{
	size_t i;
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
	Compressed_istream f(config.prot_accession2taxid);
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
	Input_stream f(config.nodesdmp);
	unsigned taxid, parent;
	while (!f.eof() && (f.getline(), !f.line.empty())) {
		if (sscanf(f.line.c_str(), "%u\t|\t%u", &taxid, &parent) != 2)
			throw std::runtime_error("Invalid nodes.dmp file format.");
		//cout << taxid << '\t' << parent << endl;
		parent_.resize(taxid + 1);
		parent_[taxid] = parent;
	}
	f.close();
}

void Taxonomy::init()
{
	task_timer timer;
	if (!config.prot_accession2taxid.empty()) {
		timer.go("Loading taxonomy");
		load();
		timer.finish();
	}

	if (!config.nodesdmp.empty()) {
		timer.go("Loading taxonomy nodes");
		load_nodes();
		timer.finish();
	}
}

void Taxonomy::get_taxids(const char *id, set<unsigned> &taxons) const
{
	const vector<string> t(tokenize(id, "\1"));
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