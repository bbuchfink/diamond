/****
DIAMOND protein aligner
Copyright (C) 2021 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include <list>
#include <vector>
#include <mutex>
#include "../sequence_set.h"
#include "../seed_histogram.h"
#include "../util/seq_file_format.h"
#include "../masking/masking.h"

struct SequenceFile;

struct SeqInfo {
	BlockId block_id;
	OId oid;
	const char* title, *qual;
	Loc len;
	Sequence source_seq, mate_seq;
};

struct Block {

	Block(Alphabet alphabet = Alphabet::STD);
	unsigned source_len(unsigned block_id) const;
	TranslatedSequence translated(size_t block_id) const;
	bool long_offsets() const;
	bool empty() const;
	void convert_to_std_alph(size_t i);
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
	Alphabet alphabet() const {
		return seqs_.alphabet();
	}
	bool fetch_seq_if_unmasked(size_t block_id, std::vector<Letter>& seq);
	void write_masked_seq(size_t block_id, const std::vector<Letter>& seq);
	DictId dict_id(size_t block, BlockId block_id, SequenceFile& db) const;
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