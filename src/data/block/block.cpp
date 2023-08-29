/****
DIAMOND protein aligner
Copyright (C) 2013-2022 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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

#include <list>
#include <atomic>
#include <thread>
#include <algorithm>
#include <iterator>
#include "../sequence_set.h"
#include "../util/seq_file_format.h"
#include "block.h"
#include "../util/sequence/sequence.h"
#include "../sequence_file.h"
#include "../basic/config.h"
#include "../dp/ungapped.h"
#define _REENTRANT
#include "../../lib/ips4o/ips4o.hpp"

using std::vector;
using std::mutex;
using std::lock_guard;
using std::array;
using std::pair;
using std::greater;
using std::string;
using std::runtime_error;
using std::numeric_limits;

Block::Block(Alphabet alphabet):
	seqs_(alphabet),
	source_seqs_(Alphabet::STD),
	unmasked_seqs_(alphabet),
	soft_masked_(false)
{
}

bool Block::empty() const {
	return seqs_.size() == 0;
}

void Block::convert_to_std_alph(size_t block_id)
{
	throw std::runtime_error("");
}

unsigned Block::source_len(unsigned block_id) const
{
	return align_mode.query_translated ? (unsigned)seqs_.reverse_translated_len(block_id * align_mode.query_contexts) : (unsigned)seqs_.length(block_id);
}

TranslatedSequence Block::translated(size_t block_id) const
{
	if (align_mode.query_translated)
		return seqs_.translated_seq(source_seqs_[block_id], block_id * align_mode.query_contexts);
	else
		return TranslatedSequence(seqs_[block_id]);
}

bool Block::long_offsets() const {
	return seqs_.raw_len() > (int64_t)std::numeric_limits<uint32_t>::max();
}

int64_t Block::push_back(const Sequence& seq, const char* id, const std::vector<char>* quals, const OId oid, const SequenceType seq_type, const int frame_mask, const bool dna_translation) {
	static const char* const OVERFLOW_ERR = "Sequences in block exceed supported maximum.";
	if (block2oid_.size() == numeric_limits<BlockId>::max())
		throw runtime_error(OVERFLOW_ERR);
	if (id)
		ids_.push_back(id, id + strlen(id));
	if (quals)
		qual_.push_back(quals->cbegin(), quals->cend());
	block2oid_.push_back(oid);
	if (seq_type == SequenceType::amino_acid || !dna_translation) {
		seqs_.push_back(seq.data(), seq.end());
		return seq.length();
	}
	else {
		if (seqs_.size() > numeric_limits<BlockId>::max() - 6)
			throw runtime_error(OVERFLOW_ERR);
		source_seqs_.push_back(seq.data(), seq.end());
		auto t = Util::Seq::translate(seq);
		const Loc min_len = config.min_orf_len((Loc)t.front().size());
		int64_t letters = 0;
		for (int j = 0; j < 6; ++j) {
			if (frame_mask & (1 << j)) {
				letters += Util::Seq::find_orfs(t[j], min_len);
				seqs_.push_back(t[j].cbegin(), t[j].cend());
			}
			else
				seqs_.fill(t[j].size(), MASK_LETTER);
		}
		return letters;
	}
}

void Block::append(const Block& b) {
	seqs_.append(b.seqs_);
	ids_.append(b.ids_);
	block2oid_.insert(block2oid_.end(), b.block2oid_.begin(), b.block2oid_.end());
}

bool Block::fetch_seq_if_unmasked(size_t block_id, std::vector<Letter>& seq) {
	if (masked_[block_id])
		return false;
	{
		lock_guard<mutex> lck(mask_lock_);
		if (masked_[block_id])
			return false;
		seq.clear();
		Sequence s = seqs_[block_id];
		std::copy(s.data(), s.end(), std::back_inserter(seq));
		return true;
	}
}

void Block::write_masked_seq(size_t block_id, const std::vector<Letter>& seq) {
	lock_guard<mutex> lck(mask_lock_);
	if (masked_[block_id])
		return;
	std::copy(seq.begin(), seq.end(), seqs_.ptr(block_id));
	masked_[block_id] = true;
}

DictId Block::dict_id(size_t block, BlockId block_id, SequenceFile& db) const
{
	string t;
	if (has_ids()) {
		const char* title = ids()[block_id];
		if (config.salltitles)
			t = title;
		else if (config.sallseqid)
			t = Util::Seq::all_seqids(title);
		else
			t = Util::Seq::seqid(title, config.short_seqids);
	}
	const Letter* seq = unmasked_seqs().empty() ? nullptr : unmasked_seqs()[block_id].data();
	double self_aln_score = 0.0;
	if (flag_any(db.flags(), SequenceFile::Flags::SELF_ALN_SCORES)) {
		if (!has_self_aln())
			throw std::runtime_error("Missing self alignment scores in Block.");
		self_aln_score = this->self_aln_score(block_id);
	}
	return db.dict_id(block, block_id, block_id2oid(block_id), seqs().length(block_id), t.c_str(), seq, self_aln_score);
}

void Block::soft_mask(const MaskingAlgo algo) {
	if (soft_masked_)
		return;
	if (soft_masking_table_.blank())
		mask_seqs(seqs_, Masking::get(), true, algo, &soft_masking_table_);
	else
		soft_masking_table_.apply(seqs_);
	soft_masked_ = true;
}

void Block::remove_soft_masking(const int template_len, const bool add_bit_mask) {
	if (!soft_masked_)
		return;
	soft_masking_table_.remove(seqs_, template_len, add_bit_mask);
	soft_masked_ = false;
}

bool Block::soft_masked() const {
	return soft_masked_;
}

size_t Block::soft_masked_letters() const {
	return soft_masking_table_.masked_letters();
}

void Block::compute_self_aln() {
	self_aln_score_.resize(seqs_.size());
	std::atomic_size_t next(0);
	auto worker = [this, &next] {
		const size_t n = this->seqs_.size();
		size_t i;
		while ((i = next++) < n) {
			this->seqs().convert_to_std_alph(i);
			this->self_aln_score_[i] = score_matrix.bitscore(self_score(this->seqs_[i]));
		}
	};
	vector<std::thread> t;
	for (int i = 0; i < config.threads_; ++i)
		t.emplace_back(worker);
	for (auto& i : t)
		i.join();
	this->seqs().alphabet() = Alphabet::STD;
}

double Block::self_aln_score(const int64_t block_id) const {
	return self_aln_score_[block_id];
}

BlockId Block::oid2block_id(OId i) const {
	if (block2oid_.back() - block2oid_.front() != seqs_.size() - 1)
		throw std::runtime_error("Block has a sparse OId range.");
	if (i < block2oid_.front() || i > block2oid_.back())
		throw std::runtime_error("OId not contained in block.");
	return BlockId(i - block2oid_.front());
}

SeqInfo Block::seq_info(const BlockId id) const {
	static const char* blank = "";
	auto mate_id = (id % 2) == 0 ? id + 1 : id - 1;
	return { id, block_id2oid(id), ids_.empty() ? nullptr : ids_[id], qual_.empty() ? blank : qual_[id],
		align_mode.query_translated ? source_seqs_.length(id) : seqs_.length(id),
		align_mode.query_translated ? source_seqs_[id] : seqs_[id],
		align_mode.query_translated && mate_id < source_seqs_.size() ? source_seqs_[mate_id] : Sequence() };
}

Block* Block::length_sorted(int threads) const {
	const BlockId n = seqs_.size();
	vector<pair<Loc, BlockId>> lengths = seqs_.lengths();
#if _MSC_FULL_VER == 191627045 || !defined(NDEBUG)
	std::sort(lengths.begin(), lengths.end(), greater<pair<Loc, BlockId>>());
#else
	ips4o::parallel::sort(lengths.begin(), lengths.end(), greater<pair<Loc, BlockId>>(), threads);
#endif
	Block *b = new Block();
	for (BlockId i = 0; i < n; ++i) {
		const BlockId j = lengths[i].second;
		b->seqs_.reserve(seqs_.length(j));
		b->ids_.reserve(ids_.length(j));
	}
	b->seqs_.finish_reserve();
	b->ids_.finish_reserve();
	b->block2oid_.reserve(n);
	for (BlockId i = 0; i < n; ++i) {
		const BlockId j = lengths[i].second;
		b->seqs_.assign(i, seqs_.ptr(j), seqs_.end(j));
		b->ids_.assign(i, ids_.ptr(j), ids_.end(j));
		b->block2oid_.push_back(block2oid_.at(j));
	}
	return b;
}

bool Block::has_ids() const {
	return !ids_.empty();
}

BlockId Block::source_seq_count() const {
	return source_seqs_.empty() ? seqs_.size() : source_seqs_.size();
}

int64_t Block::mem_size() const {
	return seqs_.mem_size() + source_seqs_.mem_size() + unmasked_seqs_.mem_size() + ids_.mem_size() + qual_.mem_size() + block2oid_.size() * sizeof(OId) + soft_masking_table_.mem_size();
}