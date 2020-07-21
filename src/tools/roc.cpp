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
#include <mutex>
#include <condition_variable>
#include <queue>
#include <thread>
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
using std::queue;
using std::thread;

typedef tuple<char, int, int, int> Family;

static map<Family, int> fam2idx;

struct FamilyMapping : public unordered_multimap<string, int> {

	FamilyMapping() {}

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

static FamilyMapping acc2fam;
static FamilyMapping acc2fam_query;
static int families;
static vector<int> fam_count;
static std::mutex mtx, mtx_out;
static std::condition_variable cdv;

struct QueryStats {
	QueryStats(const string &query, int families, const unordered_multimap<string, int>& acc2fam):
		query(query),
		count(families, 0),
		have_rev_hit(false)
	{
		if (config.output_hits) {
			query_family.insert(query_family.begin(), families, false);
			auto i = acc2fam.equal_range(query);
			for (auto j = i.first; j != i.second; ++j)
				query_family[j->second] = true;
		}
	}
	bool add(const string &sseqid, const unordered_multimap<string, int> &acc2fam) {
		if (have_rev_hit)
			return false;
		if (sseqid == last_subject)
			return false;
		if (sseqid[0] == '\\') {
			have_rev_hit = true;
			return false;
		}
		else {
			bool match_query = false;
			auto i = acc2fam.equal_range(sseqid);
			if (i.first == i.second)
				throw std::runtime_error("Accession not mapped.");
			for (auto j = i.first; j != i.second; ++j) {
				++count[j->second];
				if (config.output_hits && query_family[j->second])
					match_query = true;
			}
			last_subject = sseqid;
			return match_query;
		}
	}
	double auc1(const vector<int> &fam_count, const unordered_multimap<string, int> &acc2fam) const {
		auto i = acc2fam.equal_range(query);
		if (i.first == i.second)
			throw std::runtime_error("Query accession not mapped.");
		double r = 0.0, n = 0.0;
		for (auto j = i.first; j != i.second; ++j) {
			if (fam_count[j->second] == 0)
				continue;
			r += double(count[j->second]) / double(fam_count[j->second]);
			n += 1.0;
		}
		return n > 0.0 ? r / n : 1.0;
	}
	string query, last_subject;
	vector<int> count;
	vector<bool> query_family;
	bool have_rev_hit;
};

double query_roc(const string& buf) {
	string acc;
	Util::String::Tokenizer(buf, "\t") >> acc;
	QueryStats stats(acc, families, acc2fam_query);
	Util::String::Tokenizer tok(buf, "\t");
	while(tok.good() && !stats.have_rev_hit) {
		tok >> Util::String::Skip() >> acc;
		/*if (r.qseqid.empty() || r.sseqid.empty() || std::isnan(r.evalue) || !std::isfinite(r.evalue) || r.evalue < 0.0)
			throw std::runtime_error("Format error.");*/
		if (stats.add(acc, acc2fam) && config.output_hits)
			cout << acc << endl;
		tok.skip_to('\n');
	}
	const double a = stats.auc1(fam_count, acc2fam_query);
	if (!config.output_hits) {
		std::lock_guard<std::mutex> lock(mtx_out);
		cout << stats.query << '\t' << a << endl;
	}
	return a;
}

static queue<string*> buffers;
static bool finished = false;

static void worker() {
	string* buf;
	while (true) {
		{
			std::unique_lock<std::mutex> lock(mtx);
			while (buffers.empty() && !finished) cdv.wait(lock);
			if (buffers.empty() && finished)
				return;
			buf = buffers.front();
			buffers.pop();
		}
		query_roc(*buf);
		delete buf;
	}
}

void roc() {
	task_timer timer("Loading family mapping");
	acc2fam = FamilyMapping(config.family_map);
	timer.finish();
	message_stream << "#Mappings: " << acc2fam.size() << endl;
	message_stream << "#Families: " << fam2idx.size() << endl;

	timer.go("Loading query family mapping");
	acc2fam_query = FamilyMapping(config.family_map_query);
	timer.finish();
	message_stream << "#Mappings: " << acc2fam_query.size() << endl;
	message_stream << "#Families: " << fam2idx.size() << endl;

	families = (int)fam2idx.size();
	fam_count.insert(fam_count.begin(), families, 0);
	for (auto i = acc2fam.begin(); i != acc2fam.end(); ++i)
		++fam_count[i->second];

	vector<thread> threads;
	for (unsigned i = 0; i < std::min(config.threads_, 6u); ++i)
		threads.emplace_back(worker);

	timer.go("Processing alignments");
	TextInputFile in(config.query_file);
	string query, acc;
	size_t n = 0, queries = 0, buf_size = 0;
	string* buf = new string;
	
	while (in.getline(), !in.eof()) {
		if (in.line.empty())
			break;
		Util::String::Tokenizer(in.line, "\t") >> acc;
		if (acc != query) {
			if (!query.empty()) {
				std::lock_guard<std::mutex> lock(mtx);
				//buf_size = std::max(buf_size, buf.size());
				buffers.push(buf);
				buf = new string;
				//buf.reserve(buf_size);
				cdv.notify_one();
			}
			query = acc;
			++queries;
			if (queries % 10000 == 0)
				message_stream << "#Queries = " << queries << endl;
		}
		buf->append(in.line);
		buf->append("\n");
		++n;
	}
	{
		std::lock_guard<std::mutex> lock(mtx);
		buffers.push(buf);
		finished = true;
		cdv.notify_all();
	}

	for (thread& t : threads)
		t.join();
	
	in.close();
	timer.finish();
	message_stream << "#Records: " << n << endl;
	message_stream << "#Queries: " << queries << endl;
}