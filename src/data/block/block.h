/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

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
	StringSet& ids() {
		if (ids_.empty())
			throw std::runtime_error("Block::ids()");
		return ids_;
	}
	void reserve_ids(uint64_t letters, BlockId entries) {
		ids_.reserve(entries, letters);
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
	void append(const Block& b, bool remove_padding = false);
	SeqInfo seq_info(const BlockId id) const;
	Block* length_sorted(int threads) const;
	bool has_ids() const;
	BlockId source_seq_count() const;
	int64_t mem_size() const;
	uint64_t raw_bytes() const noexcept {
		return raw_bytes_;
	}
	template<typename T>
	void load_stats(T& stream, double microseconds) const {
		if (raw_bytes() > 0)
			stream << "Loaded " << raw_bytes() << " bytes from disk at " << ((double)raw_bytes() / MEGABYTES / microseconds * 1e6) << " MB/s" << std::endl;
	}
	BlockId oid_count() const {
		return (BlockId)block2oid_.size();
	}
	void offset_oids(OId offset);

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
	uint64_t raw_bytes_;

	friend struct SequenceFile;

};