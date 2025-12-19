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
#include <vector>
#include <mutex>
#include "../sequence_set.h"
#include "../seed_histogram.h"
#include "masking/masking.h"
#include "output/output_format.h"

struct SequenceFile;

struct Block {

	Block();
	unsigned source_len(unsigned block_id) const;
	TranslatedSequence translated(size_t block_id) const;
	bool long_offsets() const;
	bool empty() const;
	SequenceSet& seqs() {
		return seqs_;
	}
	const SequenceSet& seqs() const {
		return seqs_;
	}
	const StringSet& ids() const {
		if (ids_.empty())
			throw std::runtime_error("Block::ids()");
		return ids_;
	}
	const SequenceSet& source_seqs() const {
		return source_seqs_;
	}
	SequenceSet& unmasked_seqs() {
		return unmasked_seqs_;
	}
	const SequenceSet& unmasked_seqs() const {
		return unmasked_seqs_;
	}
	const StringSet& qual() const {
		return qual_;
	}
	SeedHistogram& hst() {
		return hst_;
	}
	OId block_id2oid(BlockId i) const {
		return block2oid_[i];
	}
	OId oid_begin() const {
		return *std::min_element(block2oid_.begin(), block2oid_.end());
		//return block2oid_.front();
	}
	OId oid_end() const {
		return *std::max_element(block2oid_.begin(), block2oid_.end()) + 1;
		//return block2oid_.back() + 1;
	}
	BlockId oid2block_id(OId i) const;
	bool fetch_seq_if_unmasked(size_t block_id, std::vector<Letter>& seq);
	void write_masked_seq(size_t block_id, const std::vector<Letter>& seq);
	DictId dict_id(size_t block, BlockId block_id, SequenceFile& db, const OutputFormat& format) const;
	void soft_mask(const MaskingAlgo algo);
	void remove_soft_masking(const int template_len, const bool add_bit_mask);
	bool soft_masked() const;
	size_t soft_masked_letters() const;
	void compute_self_aln();
	double self_aln_score(const int64_t block_id) const;
	bool has_self_aln() const {
		return (BlockId)self_aln_score_.size() == seqs_.size();
	}
	int64_t push_back(const Sequence& seq, const char* id, const std::vector<char>* quals, const OId oid, const SequenceType seq_type, const int frame_mask, const bool dna_translation = true);
	void append(const Block& b);
	SeqInfo seq_info(const BlockId id) const;
	Block* length_sorted(int threads) const;
	bool has_ids() const;
	BlockId source_seq_count() const;
	int64_t mem_size() const;

private:

	SequenceSet seqs_, source_seqs_, unmasked_seqs_;
	StringSet ids_;
	StringSet qual_;
	SeedHistogram hst_;
	std::vector<OId> block2oid_;
	std::vector<bool> masked_;
	std::vector<double> self_aln_score_;
	std::mutex mask_lock_;
	MaskingTable soft_masking_table_;
	bool soft_masked_;

	friend struct SequenceFile;

};