/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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

#include <map>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <math.h>
#include "../util/io/text_input_file.h"
#include "../basic/config.h"
#include "../util/string/tokenizer.h"
#include "../util/log_stream.h"
#include "tsv_record.h"

using std::string;
using std::cout;
using std::endl;
using std::tuple;
using std::map;
using std::unordered_multimap;
using std::vector;

struct QueryStats {
	QueryStats(const string &query, int families):
		query(query),
		count(families, 0),
		have_rev_hit(false)
	{}
	void add(const TSVRecord &r, const unordered_multimap<string, int> &acc2fam) {
		if (have_rev_hit)
			return;
		if (r.sseqid == last_subject)
			return;
		if (r.sseqid[0] == '\\')
			have_rev_hit = true;
		else {
			auto i = acc2fam.equal_range(r.sseqid);
			if (i.first == i.second)
				throw std::runtime_error("Accession not mapped.");
			for (auto j = i.first; j != i.second; ++j)
				++count[j->second];
			last_subject = r.sseqid;
		}
	}
	double auc1(const vector<int> &fam_count, const unordered_multimap<string, int> &acc2fam) const {
		auto i = acc2fam.equal_range(query);
		if (i.first == i.second)
			return -1.0;
			//throw std::runtime_error("Accession not mapped.");
		double r = 0.0, n = 0.0;
		for (auto j = i.first; j != i.second; ++j) {
			r += double(count[j->second]) / double(fam_count[j->second]);
			n += 1.0;
		}
		return r / n;
	}
	string query, last_subject;
	vector<int> count;
	bool have_rev_hit;
};

void roc() {
	typedef tuple<char, int, int, int> Family;

	TextInputFile mapping_file(config.family_map);
	string acc, domain_class;
	int fold, superfam, fam;
	map<Family, int> fam2idx;
	unordered_multimap<string, int> acc2fam;

	task_timer timer("Loading family mapping");
	while (mapping_file.getline(), !mapping_file.eof()) {
		Util::String::Tokenizer(mapping_file.line, "\t") >> Util::String::Skip() >> acc >> Util::String::Skip() >> domain_class >> fold >> superfam >> fam;
		if (acc.empty() || domain_class.empty() || domain_class.length() > 1)
			throw std::runtime_error("Format error.");

		auto i = fam2idx.insert({ Family(domain_class[0], fold, superfam, fam), (int)fam2idx.size() });
		acc2fam.insert({ acc, i.first->second });
		//cout << acc << '\t' << domain_class << '\t' << fold << '\t' << superfam << '\t' << fam << endl;
	}
	mapping_file.close();
	timer.finish();
	message_stream << "#Mappings: " << acc2fam.size() << endl;
	message_stream << "#Families: " << fam2idx.size() << endl;

	const int families = (int)fam2idx.size();
	vector<int> fam_count(families);
	for (auto i = acc2fam.begin(); i != acc2fam.end(); ++i)
		++fam_count[i->second];

	timer.go("Processing alignments");
	TextInputFile in(config.query_file);
	TSVRecord r;
	size_t n = 0, queries = 0;
	double auc1 = 0.0;
	QueryStats stats("", families);
	while (in >> r, !in.eof()) {
		if (r.qseqid != stats.query) {
			if (!stats.query.empty()) {
				const double a = stats.auc1(fam_count, acc2fam);
				if (a != -1.0) {
					auc1 += a;
					++queries;
				}
			}
			stats = QueryStats(r.qseqid, families);
		}
		if (r.qseqid.empty() || r.sseqid.empty() || std::isnan(r.evalue) || !std::isfinite(r.evalue) || r.evalue > 100.0 || r.evalue < 0.0)
			throw std::runtime_error("Format error.");
		stats.add(r, acc2fam);
		++n;
	}
	const double a = stats.auc1(fam_count, acc2fam);
	if (a != -1.0) {
		auc1 += a;
		++queries;
	}
	in.close();
	timer.finish();
	message_stream << "#Records: " << n << endl;
	message_stream << "#Queries: " << queries << endl;
	message_stream << "AUC1: " << auc1 / queries << endl;
}