/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

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
#include "data/sequence_file.h"
#include "external.h"

struct ChunkSeqs {

	ChunkSeqs(const std::string& chunk_path) :
		seq_file_(chunk_path + "bucket.tsv"),
		oid_count_(0),
		letter_count_(0)
	{
		std::atomic<int> next(0);
		seq_blocks_.resize(seq_file_.size());
		oid2seq_.resize(seq_file_.size());
		oid_range_.resize(seq_file_.size());
		//std::mutex mtx;
		auto worker = [&]() {
			int i;
			while (i = next.fetch_add(1, std::memory_order_relaxed), i < (int)seq_file_.size()) {
				std::unique_ptr<SequenceFile> in(SequenceFile::auto_create({ seq_file_[i].path }));
				in->flags() |= SequenceFile::Flags::SEQS | SequenceFile::Flags::TITLES;
				seq_blocks_[i] = in->load_seqs(INT64_MAX);
				in->close();
				Block& b = *seq_blocks_[i];
				const StringSet& ids = b.ids();
				const SequenceSet& seqs = b.seqs();
				const BlockId n = ids.size();
				std::unordered_map<int64_t, Sequence>& m = *(oid2seq_[i] = new std::unordered_map<int64_t, Sequence>());
				m.reserve(n);
				int64_t oid_min = INT64_MAX, oid_max = INT64_MIN;
				for (BlockId j = 0; j < n; ++j) {
					const int64_t oid = std::atoll(ids[j]);
					oid_min = std::min(oid_min, oid);
					oid_max = std::max(oid_max, oid);
					m[oid] = seqs[j];
				}
				{
					//std::lock_guard<std::mutex> lock(mtx);
					//range2block_.insert(std::make_pair(std::pair<int64_t, int64_t>(oid_min, oid_max), i));
					oid_range_[i] = { oid_min, oid_max };
				}
			}
			};
		std::vector<std::thread> workers;
		for (int i = 0; i < std::min(config.threads_, int(seq_file_.size())); ++i)
			workers.emplace_back(worker);
		for (auto& t : workers)
			t.join();
		int64_t oid_count = 0, letter_count = 0;
		for (const Block* b : seq_blocks_) {
			oid_count_ += b->seqs().size();
			letter_count_ += b->seqs().letters();
		}
	}

	~ChunkSeqs() {
		std::vector<std::thread> workers;
		std::atomic<int> next(0);
		auto worker = [&]() {
			int i;
			while (i = next.fetch_add(1, std::memory_order_relaxed), i < (int)seq_blocks_.size()) {
				delete seq_blocks_[i];
				delete oid2seq_[i];
			}
			};
		for (int i = 0; i < std::min(config.threads_, int(seq_blocks_.size())); ++i)
			workers.emplace_back(worker);
		for (auto& t : workers)
			t.join();
		seq_file_.remove();
	}

	size_t oids() const {
		return oid_count_;
	}

	size_t letters() const {
		return letter_count_;
	}

	size_t volumes() const {
		return seq_blocks_.size();
	}

	Sequence operator[](int64_t oid) const {
		//const auto r = range2block_.equal_range(std::make_pair(oid, oid));
		//for (auto i = r.first; i != r.second; ++i) {
		for (int i = 0; i < (int)oid_range_.size(); ++i) {
			if (oid < oid_range_[i].first || oid > oid_range_[i].second)
				continue;
			std::unordered_map<int64_t, Sequence>::const_iterator it;
			if ((it = oid2seq_[i]->find(oid)) != oid2seq_[i]->end())
				return it->second;
		}
		throw std::out_of_range("ChunkSeqs");
	}

private:

	/*struct CmpRange {
		bool operator()(const std::pair<int64_t, int64_t>& a, const std::pair<int64_t, int64_t>& b) const {
			return a.second < b.first;
		}
	};*/

	VolumedFile seq_file_;
	int64_t oid_count_, letter_count_;
	std::vector<Block*> seq_blocks_;
	//std::multimap<std::pair<int64_t, int64_t>, int, ChunkSeqs::CmpRange> range2block_;
	std::vector<std::pair<int64_t, int64_t>> oid_range_;
	std::vector<std::unordered_map<int64_t, Sequence>*> oid2seq_;

};