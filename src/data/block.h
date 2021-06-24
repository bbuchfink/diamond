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
#include "sequence_set.h"
#include "seed_histogram.h"
#include "../util/seq_file_format.h"

struct SequenceFile;

struct Block {

	Block(Alphabet alphabet);
	Block(std::list<TextInputFile>::iterator file_begin,
		std::list<TextInputFile>::iterator file_end,
		const Sequence_file_format& format,
		size_t max_letters,
		const Value_traits& value_traits,
		bool with_quals,
		bool lazy_masking = false,
		size_t modulo = 1);
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
	Partitioned_histogram& hst() {
		return hst_;
	}
	size_t block_id2oid(size_t i) const {
		return block2oid_[i];
	}
	Alphabet alphabet() const {
		return seqs_.alphabet();
	}
	bool has_ids() const {
		return !ids_.empty();
	}
	bool fetch_seq_if_unmasked(size_t block_id, std::vector<Letter>& seq);
	void write_masked_seq(size_t block_id, const std::vector<Letter>& seq);
	uint32_t dict_id(size_t block, size_t block_id, SequenceFile& db) const;

private:

	SequenceSet seqs_, source_seqs_, unmasked_seqs_;
	StringSet ids_;
	StringSet qual_;
	Partitioned_histogram hst_;
	std::vector<uint32_t> block2oid_;
	std::vector<bool> masked_;
	std::mutex mask_lock_;

	friend struct SequenceFile;

};