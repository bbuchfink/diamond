/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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