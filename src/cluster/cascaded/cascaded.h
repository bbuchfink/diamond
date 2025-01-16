/****
DIAMOND protein aligner
Copyright (C) 2016-2023 Max Planck Society for the Advancement of Science e.V.
						Benjamin Buchfink

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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
#include "../cluster.h"

namespace Cluster { 

struct Cascaded : public ClusteringAlgorithm {
	~Cascaded(){};
	void run();
	std::string get_description();
	static std::string get_key(){
		return "cascaded";
	}
};

std::vector<SuperBlockId> cascaded(std::shared_ptr<SequenceFile>& db, bool linear);
std::vector<std::string> cluster_steps(double approx_id, bool linear);
std::vector<std::string> default_round_approx_id(int steps);
std::vector<std::string> default_round_cov(int steps);
int round_ccd(int round, int round_count);

struct Callback : public Consumer {
	using Edge = Util::Algo::Edge<SuperBlockId>;
	Callback() :
		count(0)
	{}
	virtual void consume(const char* ptr, size_t n) override = 0;
	TempFile edge_file;
	int64_t count;
};

struct CallbackUnidirectional : public Callback {
	virtual void consume(const char* ptr, size_t n) override {
		const char* end = ptr + n;
		while (ptr < end) {
			const auto edge = *(Output::Format::Edge::Data*)ptr;
			ptr += sizeof(Output::Format::Edge::Data);
			if (edge.qcovhsp >= config.member_cover) {
				edge_file.write(Edge((SuperBlockId)edge.target, (SuperBlockId)edge.query, edge.evalue));
				++count;
			}
			if (edge.scovhsp >= config.member_cover) {
				edge_file.write(Edge((SuperBlockId)edge.query, (SuperBlockId)edge.target, edge.evalue));
				++count;
			}
		}
	}
};

struct CallbackBidirectional : public Callback {
	virtual void consume(const char* ptr, size_t n) override {
		const char* end = ptr + n;
		while (ptr < end) {
			const auto edge = *(Output::Format::Edge::Data*)ptr;
			ptr += sizeof(Output::Format::Edge::Data);
			if (edge.query != edge.target) {
				edge_file.write(Edge((SuperBlockId)edge.target, (SuperBlockId)edge.query, edge.evalue));
				edge_file.write(Edge((SuperBlockId)edge.query, (SuperBlockId)edge.target, edge.evalue));
				count += 2;
			}
		}
	}
};

}