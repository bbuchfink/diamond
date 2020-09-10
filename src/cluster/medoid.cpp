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

#include <limits.h>
#include <map>
#include <vector>
#include <algorithm>
#include <locale>
#include "../basic/config.h"
#include "../util/string/tokenizer.h"
#include "../util/log_stream.h"
#include "../data/reference.h"
#include "../run/workflow.h"
#include "../basic/statistics.h"
#include "../util/sequence/sequence.h"

using std::string;
using std::endl;
using std::map;
using std::pair;
using std::vector;

struct ClusterDist : public Consumer {
	virtual void consume(const char *ptr, size_t n) override {
		int query, subject, count, score;
		//double evalue;
		const char *end = ptr + n;
		while (ptr < end) {
			//if (sscanf(ptr, "%i\t%i\n%n", &query, &subject, &count) != 2)
			if (sscanf(ptr, "%i\t%i\t%i\n%n", &query, &subject, &score, &count) != 3)
				throw runtime_error("Cluster format error.");
			ptr += count;
			//std::cout << query << '\t' << subject << '\t' << evalue << endl;
			if (query != subject) {
				sum[query] += score;
				++counts[query];
			}
		}
	}
	map<int, uint64_t> sum;
	map<int, int> counts;
};

size_t get_medoid(DatabaseFile *db, const BitVector &filter, size_t n, Sequence_set *seqs) {
	statistics.reset();
	config.command = Config::blastp;
	config.no_self_hits = true;
	config.output_format = { "6", "qnum", "snum", "score" };
	config.swipe_all = true;
	config.max_evalue = 100.0;
	config.freq_sd = 0;
	config.max_alignments = SIZE_MAX;
	config.algo = 0;
	//config.ext = Config::swipe;
	score_matrix.set_db_letters(1);

	Workflow::Search::Options opt;
	opt.db = db;
	opt.self = true;
	ClusterDist d;
	opt.consumer = &d;
	opt.db_filter = &filter;

	Workflow::Search::run(opt);
	
	int max_i = -1;
	uint64_t max_s = 0;
	for (pair<int, uint64_t> i : d.sum) {
		//std::cout << i.first << '\t' << i.second  << '\t' << (*seqs)[i.first].length() << endl;
		//const double sum = i.second + (n - 1 - d.counts[i.first])*100.0;
		if (i.second > max_s) {
			max_s = i.second;
			max_i = i.first;
		}
	}
	//std::cout << "=============" << endl;
	//return max_i == -1 ? std::find(filter.begin(), filter.end(), true) - filter.begin() : max_i;
	return 0;
}

int get_acc2idx(const string& acc, const map<string, size_t>& acc2idx) {
	static std::locale loc;
	if (std::isdigit(acc[0], loc))
		return atoi(acc.c_str());
	else {
		auto i = acc2idx.find(acc);
		return i == acc2idx.end() ? -1 : (int)i->second;
	}
}

void get_medoids_from_tree() {
	const size_t CLUSTER_COUNT = 1000;
	config.masking = false;
	DatabaseFile *db = DatabaseFile::auto_create_from_fasta();
	message_stream << "#Sequences: " << db->ref_header.sequences << endl;
	size_t n = db->ref_header.sequences;

	Sequence_set *seqs;
	String_set<char, '\0'> *ids;
	db->load_seqs(nullptr, SIZE_MAX, &seqs, &ids, true);

	map<int, int> parent;
	map<string, int> acc2idx;
	for (size_t i = 0; i < n; ++i) {
		const string id = (*ids)[i];
		parent[i] = i;
		acc2idx[id] = i;
		acc2idx[std::to_string(i)] = i;
	}
	
	TextInputFile tree_in(config.tree_file);
	int p;
	string c1, c2;
	while (tree_in.getline(), !tree_in.eof() && n > CLUSTER_COUNT) {
		Util::String::Tokenizer(tree_in.line, "\t") >> p >> c1 >> c2;
		parent[acc2idx.find(c1) == acc2idx.end() ? atoi(c1.c_str()) : acc2idx[c1]] = p;
		parent[acc2idx.find(c2) == acc2idx.end() ? atoi(c2.c_str()) : acc2idx[c2]] = p;
		parent[p] = p;
		--n;
	}

	n = db->ref_header.sequences;
	map<int, vector<int>> clusters;
	for (pair<int, int> i : parent) {
		while (parent[parent[i.first]] != parent[i.first])
			parent[i.first] = parent[parent[i.first]];
		if(i.first < (int)n)
			clusters[parent[i.first]].push_back(i.first);
	}

	BitVector filter(db->ref_header.sequences);
	OutputFile out(config.output_file);
	for (const pair<const int, vector<int>> &i : clusters) {
		/*for (const string &acc : i.second)
			std::cout << acc << ' ';
		std::cout << endl;*/
		size_t medoid;
		if (i.second.size() == 1)
			medoid = i.second.front();
		else {
			filter.reset();
			for (const int& acc : i.second)
				filter.set(acc);
			medoid = get_medoid(db, filter, i.second.size(), seqs);
		}
		const string id = string((*ids)[medoid]) + ' ' + std::to_string(i.second.size());
		Util::Sequence::format((*seqs)[medoid], id.c_str(), nullptr, out, "fasta", amino_acid_traits);
	}
	out.close();

	db->close_and_delete();
	delete db;
	delete seqs;
	delete ids;
}