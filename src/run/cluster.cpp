/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include <algorithm>
#include <stdexcept>
#include <string>
#include <stdio.h>
#include <memory>
#include <fstream>
#include "../util/system/system.h"
#include "../util/util.h"
#include "../basic/config.h"
#include "../data/reference.h"
#include "workflow.h"
#include "../util/io/consumer.h"
#include "../util/algo/algo.h"
#include "../basic/statistics.h"
#include "../util/log_stream.h"

using namespace std;

namespace Workflow { namespace Cluster {

struct Neighbors : public vector<vector<int>>, public Consumer {
	Neighbors(size_t n):
		vector<vector<int>>(n)
	{}
	virtual void consume(const char *ptr, size_t n) override {
		int query, subject, count;
		float qcov, scov;
		const char *end = ptr + n;
		while (ptr < end) {
			//if (sscanf(ptr, "%i\t%i\n%n", &query, &subject, &count) != 2)
			if (sscanf(ptr, "%i\t%i\t%f\t%f\n%n", &query, &subject, &qcov, &scov, &count) != 4)
				throw runtime_error("Cluster format error.");
			ptr += count;
			//cout << query << '\t' << subject << '\t' << qcov << '\t' << scov << '\t' << endl;
			(*this)[query].push_back(subject);
		}
	}
};

vector<bool> rep_bitset(const vector<int> &centroid, const vector<bool> *superset = nullptr) {
	vector<bool> r(centroid.size());
	for (int c : centroid)
		if(!superset || (*superset)[c])
			r[c] = true;
	return r;
}

vector<int> cluster(DatabaseFile &db, const vector<bool> *filter) {
	statistics.reset();
	config.command = Config::blastp;
	config.no_self_hits = true;
	//config.output_format = { "6", "qnum", "snum" };
	config.output_format = { "6", "qnum", "snum", "qcovhsp", "scovhsp" };
	config.query_cover = 80;
	config.subject_cover = 80;
	config.algo = 0;
	config.index_mode = 0;
	config.freq_sd = 0;

	Workflow::Search::Options opt;
	opt.db = &db;
	opt.self = true;
	Neighbors nb(db.ref_header.sequences);
	opt.consumer = &nb;
	opt.db_filter = filter;

	Workflow::Search::run(opt);

	return greedy_vortex_cover(nb);
}

void run() {
	if (config.database == "")
		throw runtime_error("Missing parameter: database file (--db/-d)");
	config.command = Config::makedb;
	unique_ptr<DatabaseFile> db(DatabaseFile::auto_create_from_fasta());

	vector<int> centroid1(cluster(*db, nullptr));
	vector<bool> rep1(rep_bitset(centroid1));
	size_t n_rep1 = count_if(rep1.begin(), rep1.end(), [](bool x) { return x; });
	message_stream << "Clustering step 1 complete. #Input sequences: " << centroid1.size() << " #Clusters: " << n_rep1 << endl;

	config.mode_more_sensitive = true;
	vector<int> centroid2(cluster(*db, &rep1));
	vector<bool> rep2(rep_bitset(centroid2, &rep1));
	for (size_t i = 0; i < centroid2.size(); ++i)
		if (!rep1[i])
			centroid2[i] = centroid2[centroid1[i]];
	message_stream << "Clustering step 2 complete. #Input sequences: " << n_rep1 << " #Clusters: " << count_if(rep2.begin(), rep2.end(), [](bool x) { return x; }) << endl;

	task_timer timer("Generating output");
	ostream *out = config.output_file.empty() ? &cout : new ofstream(config.output_file.c_str());
	vector<char> seq;
	string id;
	db->seek_direct();
	for (size_t i = 0; i < db->ref_header.sequences; ++i) {
		db->read_seq(id, seq);
		(*out) << blast_id(id) << endl;
	}

	db->close();
	if (out != &cout) delete out;
}

}}