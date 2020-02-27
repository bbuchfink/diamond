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
#include <limits>
#include "../util/system/system.h"
#include "../util/util.h"
#include "../basic/config.h"
#include "../data/reference.h"
#include "workflow.h"
#include "disjoint_set.h"
#include "../util/io/consumer.h"
#include "../util/algo/algo.h"
#include "../basic/statistics.h"
#include "../util/log_stream.h"
#include "../dp/dp.h"
#include "../basic/masking.h"

using namespace std;

namespace Workflow { namespace Cluster {

struct Neighbors : public vector<vector<int>>, public Consumer {
	Neighbors(size_t n):
		vector<vector<int>>(n)
	{}
	virtual void consume(const char *ptr, size_t n) override {
		int query, subject, count;
		float qcov, scov, bitscore;
		const char *end = ptr + n;
		while (ptr < end) {
			//if (sscanf(ptr, "%i\t%i\n%n", &query, &subject, &count) != 2)
			if (sscanf(ptr, "%i\t%i\t%f\t%f\t%f\n%n", &query, &subject, &qcov, &scov, &bitscore, &count) != 5)
				throw runtime_error("Cluster format error.");
			ptr += count;
			//cout << query << '\t' << subject << '\t' << qcov << '\t' << scov << '\t' << endl;
			(*this)[query].push_back(subject);
			edges.push_back({ query, subject, (int)bitscore });
		}
	}
	vector<Util::Algo::Edge> edges;
};

class NeighborStream : public Consumer {
	LazyDisjointSet<size_t>* disjointSet;
	virtual void consume(const char *ptr, size_t n) override {
		size_t query, subject, count;
		float qcov, scov, bitscore;
		const char *end = ptr + n;
		while (ptr < end) {
			//if (sscanf(ptr, "%i\t%i\n%n", &query, &subject, &count) != 2)
			if (sscanf(ptr, "%lu\t%lu\t%f\t%f\t%f\n%ln", &query, &subject, &qcov, &scov, &bitscore, &count) != 5)
				throw runtime_error("Cluster format error.");
			disjointSet->merge(query, subject);
			ptr += count;
		}
	}
public:
	NeighborStream(size_t n){
		disjointSet = new LazyDisjointIntegralSet<size_t>(n); 
	}
	~NeighborStream(){
		delete disjointSet;
	}
	vector<unordered_set<size_t>> getListOfSets(){
		return disjointSet->getListOfSets();
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
	config.output_format = { "6", "qnum", "snum", "qcovhsp", "scovhsp", "bitscore" };
	config.query_cover = 80;
	config.subject_cover = 80;
	config.algo = 0;
	config.index_mode = 0;
	config.freq_sd = 0;
	config.max_alignments = numeric_limits<uint64_t>::max();

	Workflow::Search::Options opt;
	opt.db = &db;
	opt.self = true;
	Neighbors nb(db.ref_header.sequences);
	opt.consumer = &nb;
	opt.db_filter = filter;

	Workflow::Search::run(opt);
	
	return Util::Algo::greedy_vortex_cover(nb);
}

void run_two_step_clustering(DatabaseFile& db){
	const size_t seq_count = db.ref_header.sequences;

	config.min_id = 70;
	vector<int> centroid1(cluster(db, nullptr));
	vector<bool> rep1(rep_bitset(centroid1));
	size_t n_rep1 = count_if(rep1.begin(), rep1.end(), [](bool x) { return x; });
	message_stream << "Clustering step 1 complete. #Input sequences: " << centroid1.size() << " #Clusters: " << n_rep1 << endl;

	config.mode_more_sensitive = true;
	config.min_id = 0;
	vector<int> centroid2(cluster(db, &rep1));
	vector<bool> rep2(rep_bitset(centroid2, &rep1));
	for (size_t i = 0; i < centroid2.size(); ++i)
		if (!rep1[i])
			centroid2[i] = centroid2[centroid1[i]];
	message_stream << "Clustering step 2 complete. #Input sequences: " << n_rep1 << " #Clusters: " << count_if(rep2.begin(), rep2.end(), [](bool x) { return x; }) << endl;

	task_timer timer("Generating output");
	Sequence_set *rep_seqs;
	String_set<0> *rep_ids;
	vector<unsigned> rep_database_id, rep_block_id(seq_count);
	db.rewind();
	db.load_seqs(rep_database_id, (size_t)1e11, &rep_seqs, &rep_ids, true, &rep2);
	for (size_t i = 0; i < rep_database_id.size(); ++i)
		rep_block_id[rep_database_id[i]] = (unsigned)i;

	ostream *out = config.output_file.empty() ? &cout : new ofstream(config.output_file.c_str());
	vector<char> seq;
	string id;
	db.seek_direct();
	Hsp hsp;
	size_t n;
	out->precision(3);

	for (int i = 0; i < (int)db.ref_header.sequences; ++i) {
		db.read_seq(id, seq);
		const unsigned r = rep_block_id[centroid2[i]];

		(*out) << blast_id(id) << '\t'
			<< blast_id((*rep_ids)[r].c_str()) << '\t';		
		
		if ((int)i == centroid2[i])
			(*out) << "100\t100\t100\t0" << endl;
		else {
			Masking::get().bit_to_hard_mask(seq.data(), seq.size(), n);
			smith_waterman(sequence(seq), (*rep_seqs)[r], hsp);
			(*out) << hsp.id_percent() << '\t'
				<< hsp.query_cover_percent((unsigned)seq.size()) << '\t'
				<< hsp.subject_cover_percent((unsigned)(*rep_seqs)[r].length()) << '\t'
				<< score_matrix.bitscore(hsp.score) << endl;
		}
	}

	if (out != &cout) delete out;
	delete rep_seqs;
	delete rep_ids;
}


void run_transitive_closure_clustering(DatabaseFile& db){
	// testing
	// LazyDisjointIntegralSet<uint32_t> test(5);
	// test.merge(1,3);
	// test.merge(3,4);
	// test.merge(0,2);
	// auto los = test.getListOfSets();
	// for(auto& set: los){
	// 	for(auto item:set){
	// 		cout<< " " << item;
	// 	}
	// 	cout<<endl;
	// }
	// unordered_set<string> s ( {"red","green","blue","rain", "sun", "snow"} );
	// LazyDisjointTypeSet<string> ts(&s);

	// ts.merge("red", "green");
	// ts.merge("snow", "rain");
	// ts.merge("sun", "rain");
	// auto lo = ts.getListOfSets();
	// for(auto& set: lo){
	// 	cout<< "set :";
	// 	for(auto item:set){
	// 		cout<< " " << item;
	// 	}
	// 	cout<<endl;
	// }

	const size_t seq_count = db.ref_header.sequences;
	statistics.reset();
	config.command = Config::blastp;
	config.no_self_hits = true;
	//config.output_format = { "6", "qnum", "snum" };
	config.output_format = { "6", "qnum", "snum", "qcovhsp", "scovhsp", "bitscore" };
	config.query_cover = 80;
	config.subject_cover = 80;
	config.algo = 0;
	config.index_mode = 0;
	config.freq_sd = 0;
	config.max_alignments = numeric_limits<uint64_t>::max();

	Workflow::Search::Options opt;
	opt.db = &db;
	opt.self = true;
	NeighborStream ns(db.ref_header.sequences);
	opt.consumer = &ns;
	opt.db_filter = nullptr;

	Workflow::Search::run(opt);

	auto clusterResult = ns.getListOfSets();
	cout << "Found "<<clusterResult.size()<<" clusters"<<endl;
	// for(auto& set: clusterResult){
	// 	cout<< "set :";
	// 	for(auto item : set){
	// 		cout<< " " << item;
	// 	}
	// 	cout<<endl;
	// }
}

void run() {
	if (config.database == "")
		throw runtime_error("Missing parameter: database file (--db/-d)");
	config.command = Config::makedb;
	unique_ptr<DatabaseFile> db(DatabaseFile::auto_create_from_fasta());
	switch(config.cluster_algo){
		case Config::multi_step:
			run_two_step_clustering(*db);
			break;
		case Config::transitive_closure:
			run_transitive_closure_clustering(*db);
			break;
	}
	db->close();
}
}}
