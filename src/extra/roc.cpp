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
#include <map>
#include <set>
#include "../util/binary_file.h"
#include "../basic/config.h"
#include "match_file.h"
#include "../util/seq_file_format.h"

using std::map;
using std::set;

const int roc_from = -10, roc_to = 1;
const size_t roc_steps = (roc_to - roc_from + 1) * 9;

struct Superfamily
{
	bool operator==(const Superfamily &x) const
	{
		return cl == x.cl && fold == x.fold && superfamily == x.superfamily;
	}
	char cl;
	unsigned fold, superfamily;
};

map<char, map<unsigned, map<unsigned, unsigned> > > superfamilies;
map<string, Superfamily> subjects;
set<pair<string, string> > target;
size_t n_targets = 0, tp = 0, fp = 0;
double max_ev = 0;

void query_roc(Superfamily superfamily, const match_file::mcont &matches, Numeric_vector<double> &coverage, Numeric_vector<double> &errors)
{
	coverage = Numeric_vector<double>(roc_steps);
	errors = Numeric_vector<double>(roc_steps);
	match_file::mcont::const_iterator i = matches.begin();
	size_t idx = 0;
	for (int exp = roc_from; exp <= roc_to; ++exp)
		for (int factor = 2; factor <= 10; ++factor) {
			const double ev = pow(10.0, exp)*factor;
			if (idx != 0) {
				coverage[idx] = coverage[idx - 1];
				errors[idx] = errors[idx - 1];
			}
			while (i < matches.end()) {
				if (i != matches.begin() && i->subject == (i - 1)->subject) {
					++i;
					continue;
				}
				if (i->expect > ev)
					break;
				if (subjects.find(i->subject) != subjects.end() && subjects[i->subject] == superfamily) {
					++coverage[idx];
					if (target.find(pair<string, string>(i->query, i->subject)) != target.end()) {
						max_ev = std::max(max_ev, i->expect);
						++n_targets;
					}
				}
				else {
					++errors[idx];
					++fp;
				}
				++i;
			}
			++idx;
		}
	coverage /= superfamilies[superfamily.cl][superfamily.fold][superfamily.superfamily];
}

void roc()
{
	vector<char> id;
	vector<Letter> seq;

	Input_stream seqStream(config.query_file);
	match_file file1(config.match_file1.c_str());
	match_file::mcont v1;

	Numeric_vector<double> coverage(roc_steps), errors(roc_steps), c2(roc_steps), e2(roc_steps);
	size_t queries = 0;
	FASTA_format format;
	while (format.get_seq(id, seq, seqStream)) {
		string id2(id.data(), id.size());
		++queries;

		char name[32];
		unsigned family;
		Superfamily superfamily;
		if (sscanf(id2.c_str(), "%s %c.%u.%u.%u", name, &superfamily.cl, &superfamily.fold, &superfamily.superfamily, &family) != 5)
			throw std::runtime_error("Format error");
		++superfamilies[superfamily.cl][superfamily.fold][superfamily.superfamily];
		subjects[name] = superfamily;
	}

	if (!config.match_file2.empty()) {
		Input_stream target_file(config.match_file2.c_str());
		while (target_file.getline(), !target_file.eof()) {
			char q[16], s[128];
			float b;
			if (sscanf(target_file.line.c_str(), "%s %s %f", q, s, &b) != 3)
				throw std::runtime_error("Format error");
			if (s[0] == 'd' || s[0] == 'g')
				target.insert(pair<string, string>(q, s));
		}
	}

	while (file1.get_read(v1, blast_tab_format())) {
		query_roc(subjects[v1[0].query], v1, c2, e2);
		coverage += c2;
		errors += e2;
	}
	
	coverage /= (double)queries;
	errors /= (double)queries;
	cout << queries << " Sequences." << endl;
	cout << coverage << endl;
	cout << errors << endl;
	/*for (int exp = roc_from; exp <= roc_to; ++exp)
		for (int factor = 2; factor <= 10; ++factor)
			cout << pow(10.0, exp)*factor << endl;*/

	cout << endl;
	cout << "Targets = " << n_targets << " / " << target.size() << " (" << percentage<double,size_t>(n_targets, target.size()) << "%)" << endl;
	cout << "max ev = " << max_ev << endl;
	cout << "False positives = " << fp << endl;
}