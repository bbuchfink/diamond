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
#include "../util/string/tokenizer.h"

using std::string;
using std::cout;
using std::endl;
using std::tuple;
using std::map;
using std::unordered_multimap;
using std::vector;

typedef tuple<char, int, int, int> Family;

static map<Family, int> fam2idx;

struct FamilyMapping : public unordered_multimap<string, int> {

	FamilyMapping(const string &file_name) {

		TextInputFile mapping_file(file_name);
		string acc, domain_class;
		int fold, superfam, fam;

		while (mapping_file.getline(), !mapping_file.eof()) {
			Util::String::Tokenizer(mapping_file.line, "\t") >> Util::String::Skip() >> acc >> Util::String::Skip() >> domain_class >> fold >> superfam >> fam;
			if (acc.empty() || domain_class.empty() || domain_class.length() > 1)
				throw std::runtime_error("Format error.");

			auto i = fam2idx.insert({ Family(domain_class[0], fold, superfam, fam), (int)fam2idx.size() });
			if (config.cut_bar) {
				size_t j = acc.find_last_of('|');
				if (j != string::npos)
					acc = acc.substr(j + 1);
			}
			insert({ acc, i.first->second });
			//cout << acc << '\t' << domain_class << '\t' << fold << '\t' << superfam << '\t' << fam << endl;
		}
	}

};

struct QueryStats {
	QueryStats(const string &query, int families):
		query(query),
		count(families, 0),
		have_rev_hit(false)
	{}
	void add(const string &sseqid, const unordered_multimap<string, int> &acc2fam) {
		if (have_rev_hit)
			return;
		if (sseqid == last_subject)
			return;
		if (sseqid[0] == '\\')
			have_rev_hit = true;
		else {
			auto i = acc2fam.equal_range(sseqid);
			if (i.first == i.second)
				throw std::runtime_error("Accession not mapped.");
			for (auto j = i.first; j != i.second; ++j)
				++count[j->second];
			last_subject = sseqid;
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
	task_timer timer("Loading family mapping");
	FamilyMapping acc2fam(config.family_map);
	timer.finish();
	message_stream << "#Mappings: " << acc2fam.size() << endl;
	message_stream << "#Families: " << fam2idx.size() << endl;

	timer.go("Loading query family mapping");
	FamilyMapping acc2fam_query(config.family_map_query);
	timer.finish();
	message_stream << "#Mappings: " << acc2fam_query.size() << endl;
	message_stream << "#Families: " << fam2idx.size() << endl;

	const int families = (int)fam2idx.size();
	vector<int> fam_count(families);
	for (auto i = acc2fam.begin(); i != acc2fam.end(); ++i)
		++fam_count[i->second];

	timer.go("Processing alignments");
	TextInputFile in(config.query_file);
	string qseqid, sseqid;
	size_t n = 0, queries = 0;
	double auc1 = 0.0;
	QueryStats stats("", families);
	while (in.getline(), !in.eof()) {
		if (in.line.empty())
			break;
		Util::String::Tokenizer(in.line, "\t") >> qseqid >> sseqid;
		if (qseqid != stats.query) {
			if (!stats.query.empty()) {
				const double a = stats.auc1(fam_count, acc2fam_query);
				if (a != -1.0) {
					auc1 += a;
					++queries;
					cout << stats.query << '\t' << a << endl;
				}
			}
			stats = QueryStats(qseqid, families);
		}
		/*if (r.qseqid.empty() || r.sseqid.empty() || std::isnan(r.evalue) || !std::isfinite(r.evalue) || r.evalue < 0.0)
			throw std::runtime_error("Format error.");*/
		stats.add(sseqid, acc2fam);
		++n;
	}
	const double a = stats.auc1(fam_count, acc2fam_query);
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