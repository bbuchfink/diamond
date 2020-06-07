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

#pragma once
#include <algorithm>
#include <stdexcept>
#include <stdio.h>
#include <memory>
#include <fstream>
#include <limits>
#include "../util/system/system.h"
#include "../util/util.h"
#include "../basic/config.h"
#include "../data/reference.h"
#include "../run/workflow.h"
#include "../util/io/consumer.h"
#include "../util/algo/algo.h"
#include "../basic/statistics.h"
#include "../util/log_stream.h"
#include "../dp/dp.h"
#include "../basic/masking.h"
#include "cluster.h"

using namespace std;

namespace Workflow { namespace Cluster{ 

class MultiStep : public ClusteringAlgorithm {
private:
	vector<bool> rep_bitset(const vector<int> &centroid, const vector<bool> *superset = nullptr);
	vector<int> cluster(DatabaseFile &db, const vector<bool> *filter);
public:
	void run();
	string get_key();
	string get_description();
};

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

}}
