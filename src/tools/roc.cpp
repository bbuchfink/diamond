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

#include <iostream>
#include <atomic>
#include <array>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <vector>
#include <math.h>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <thread>
#include <limits.h>
#include <float.h>
#include <sstream>
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
using std::unordered_set;
using std::vector;
using std::queue;
using std::thread;
using std::set;
using std::array;
using std::stringstream;

typedef tuple<char, int> Fold;
typedef tuple<char, int, int, int> Family;

static map<Family, int> fam2idx;
static map<int, Fold> fam2fold;

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
			fam2fold[i.first->second] = Fold(domain_class[0], fold);
			//cout << acc << '\t' << domain_class << '\t' << fold << '\t' << superfam << '\t' << fam << endl;
		}
	}

};

static FamilyMapping acc2fam;
static FamilyMapping acc2fam_query;
static int families;
static vector<int> fam_count;
static std::mutex mtx, mtx_out, mtx_hist;
static std::condition_variable cdv;
static bool get_roc;
static std::atomic<size_t> query_with_fp(0);

static double coverage(int count, int family) {
	const int n = fam_count[family];
	if (n == 0)
		return 1.0;
	return double(count) / double(n);
}

struct Histogram {

	static constexpr double MAX_EV = 10000.0;

	Histogram() :
		bin_offset(int(-std::floor((double)DBL_MIN_EXP* log(2.0)* config.log_evalue_scale))),
		bin_count(bin_offset + int(std::round(std::log(MAX_EV) * config.log_evalue_scale)) + 1),
		false_positives(bin_count, 0),
		coverage(bin_count, 0.0)
	{}

	int bin(double evalue) {
		if (evalue == 0.0)
			return 0;
		int bin = int(std::round(std::log(evalue) * config.log_evalue_scale)) + bin_offset;
		bin = std::max(bin, 0);
		if (bin < 0 || bin >= bin_count)
			throw std::runtime_error("Evalue exceeds binning range.");
		return bin;
	}

	Histogram& operator+=(const Histogram& h) {
		for (int i = 0; i < bin_count; ++i) {
			false_positives[i] += h.false_positives[i];
			coverage[i] += h.coverage[i];
		}
		return *this;
	}

	friend std::ostream& operator<<(std::ostream& os, const Histogram& h) {
		for (int i = 0; i < h.bin_count; ++i)
			os << (h.coverage[i] / config.query_count) << '\t' << ((double)h.false_positives[i] / config.query_count) << endl;
		return os;
	}

	int bin_offset, bin_count;
	vector<size_t> false_positives;
	vector<double> coverage;

};

static Histogram histogram;

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
		if (!config.no_forward_fp) {
			auto i = acc2fam.equal_range(query);
			for (auto j = i.first; j != i.second; ++j)
				query_fold.insert(fam2fold[j->second]);
		}
		if (get_roc) {
			auto i = acc2fam.equal_range(query);
			int n = 0;
			for (auto j = i.first; j != i.second; ++j) {
				const int next = (int)family_idx.size();
				family_idx[j->second] = next;
				++n;
			}
			false_positives.insert(false_positives.end(), histogram.bin_count, 0);
			true_positives.reserve(n);
			for (int i = 0; i < n; ++i)
				true_positives.emplace_back(histogram.bin_count, 0);
		}
	}

	void add_family_hit(int family, double evalue) {
		if (!get_roc)
			return;
		auto it = family_idx.find(family);
		if (it == family_idx.end())
			return;
		++true_positives[it->second][histogram.bin(evalue)];
	}

	enum { UNKNOWN = 0, TP = 1, FP = 2 };

	int add(const string &sseqid, double evalue, const unordered_multimap<string, int> &acc2fam) {
		if (have_rev_hit && !get_roc)
			return UNKNOWN;
		if (sseqid == last_subject)
			return UNKNOWN;
		if (config.check_multi_target) {
			if (previous_targets.find(sseqid) != previous_targets.end())
				return false;
			previous_targets.insert(sseqid);
		}
		last_subject = sseqid;
		if (sseqid[0] == '\\') {
			have_rev_hit = true;
			if (get_roc)
				++false_positives[histogram.bin(evalue)];
			return FP;
		}
		else {
			bool match_query = false;
			auto i = acc2fam.equal_range(sseqid);
			if (i.first == i.second)
				throw std::runtime_error("Accession not mapped.");
			bool same_fold = false;

			for (auto j = i.first; j != i.second; ++j) {
				if (!have_rev_hit)
					++count[j->second];
				add_family_hit(j->second, evalue);
				if (config.output_hits && query_family[j->second])
					match_query = true;
				if (!config.no_forward_fp && query_fold.find(fam2fold[j->second]) != query_fold.end())
					same_fold = true;
			}
			if (!config.no_forward_fp && !same_fold) {
				have_rev_hit = true;
				if(get_roc)
					++false_positives[histogram.bin(evalue)];
				return FP;
			}
			return match_query ? TP : UNKNOWN;
		}
	}

	double auc1(const vector<int> &fam_count, const unordered_multimap<string, int> &acc2fam) const {
		auto i = acc2fam.equal_range(query);
		if (i.first == i.second)
			throw std::runtime_error("Query accession not mapped.");
		double r = 0.0, n = 0.0;
		for (auto j = i.first; j != i.second; ++j) {
			r += coverage(count[j->second], j->second);
			n += 1.0;
		}
		return r / n;
	}

	int family_count() const {
		return (int)family_idx.size();
	}

	void update_hist(Histogram& hist, const vector<int>& fam_count) {
		int t = 0;
		for (int i = 0; i < histogram.bin_count; ++i) {
			t += false_positives[i];
			hist.false_positives[i] += (size_t)t;
		}
		const int n = family_count();
		vector<int> s(n, 0);
		for (int i = 0; i < histogram.bin_count; ++i) {
			double cov = 0.0;
			for(auto it = family_idx.begin(); it != family_idx.end(); ++it) {
				s[it->second] += true_positives[it->second][i];
				cov += coverage(s[it->second], it->first);
			}
			hist.coverage[i] += cov / n;
		}
	}

	string query, last_subject;
	vector<int> count;
	vector<bool> query_family;
	set<Fold> query_fold;
	map<int, int> family_idx;
	vector<int> false_positives;
	unordered_set<string> previous_targets;
	vector<vector<int>> true_positives;
	bool have_rev_hit;

};

double query_roc(const string& buf, Histogram& hist) {
	string query, acc, line;
	Util::String::Tokenizer(buf, "\t") >> query;
	QueryStats stats(query, families, acc2fam_query);
	Util::String::Tokenizer tok(buf, "\n");
	double evalue = 0.0;
	while(tok.good() && (!stats.have_rev_hit || get_roc)) {
		tok >> line;
		if (line.empty())
			break;
		Util::String::Tokenizer tok2(line, "\t");
		tok2 >> Util::String::Skip() >> acc;
		if (get_roc)
			tok2 >> evalue;
		/*if (r.qseqid.empty() || r.sseqid.empty() || std::isnan(r.evalue) || !std::isfinite(r.evalue) || r.evalue < 0.0)
			throw std::runtime_error("Format error.");*/
		int c = stats.add(acc, evalue, acc2fam);
		if ((c == QueryStats::TP && config.output_hits) || (c == QueryStats::FP && config.output_fp))
			cout << line << endl;

	}
	const double a = stats.auc1(fam_count, acc2fam_query);
	if (get_roc)
		stats.update_hist(hist, fam_count);
	if (!config.output_hits && !config.output_fp) {
		std::lock_guard<std::mutex> lock(mtx_out);
		cout << stats.query << '\t' << a << endl;
	}
	if (stats.have_rev_hit)
		++query_with_fp;
	return a;
}

static queue<string*> buffers;
static bool finished = false;

static void worker() {
	string* buf;
	Histogram hist;
	while (true) {
		{
			std::unique_lock<std::mutex> lock(mtx);
			while (buffers.empty() && !finished) cdv.wait(lock);
			if (buffers.empty() && finished) {
				std::lock_guard<std::mutex> lock(mtx_hist);
				histogram += hist;
				return;
			}
			buf = buffers.front();
			buffers.pop();
		}
		query_roc(*buf, hist);
		delete buf;
	}
	{
		std::lock_guard<std::mutex> lock(mtx_hist);
		histogram += hist;
	}
}

void roc() {
	histogram = Histogram();
	get_roc = !config.roc_file.empty();

	if (config.family_map.empty())
		throw std::runtime_error("Missing option: --family-map");

	if (config.family_map_query.empty())
		throw std::runtime_error("Missing option: --family-map-query");

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
		if (config.family_cap == 0)
			++fam_count[i->second];
		else
			fam_count[i->second] = config.family_cap;

	vector<thread> threads;
	for (int i = 0; i < std::min(config.threads_, 6); ++i)
		threads.emplace_back(worker);

	timer.go("Processing alignments");
	TextInputFile in(config.single_query_file());
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
			while (buffers.size() > 100)
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
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

	if (get_roc) {
		std::ofstream out(config.roc_file);
		out << histogram;
	}

	message_stream << "#Records: " << n << endl;
	message_stream << "#Queries: " << queries << endl;
	message_stream << "#Queries w/ FP: " << query_with_fp << endl;
}