/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include "../cluster.h"

extern const double CASCADED_ROUND_MAX_EVALUE;

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
bool is_linclust(const std::vector<std::string>& steps);
std::vector<std::string> default_round_approx_id(int steps);
std::vector<std::string> default_round_cov(int steps);
int round_ccd(int round, int round_count, bool linear);

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