/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
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

#pragma once
#include <memory>
#include <unordered_set>
#include "../basic/value.h"
#include "../stats/score_matrix.h"
#include "../basic/sequence.h"
#include "../lib/blast/blast_seg.h"
#include "../data/string_set.h"
#include "../util/enum.h"
#include "def.h"

struct MaskingTable;
struct SequenceSet;

struct Masking
{
	Masking(const ScoreMatrix&score_matrix);
	~Masking();
	size_t operator()(Letter *seq, size_t len, const MaskingAlgo algo, const size_t block_id, MaskingTable* table = nullptr) const;
	void mask_bit(Letter *seq, size_t len) const;
	void bit_to_hard_mask(Letter *seq, size_t len, size_t &n) const;
	void remove_bit_mask(Letter *seq, size_t len) const;
	static const Masking& get()
	{
		return *instance;
	}
	static std::unique_ptr<Masking> instance;
	static const int8_t bit_mask;
private:
	enum { size = 64 };
	float likelihoodRatioMatrixf_[size][size], *probMatrixPointersf_[size];
	Letter mask_table_x_[size], mask_table_bit_[size];
	SegParameters* blast_seg_;
};

size_t mask_seqs(SequenceSet &seqs, const Masking &masking, bool hard_mask, const MaskingAlgo algo, MaskingTable* table = nullptr);

template<>
struct EnumTraits<MaskingAlgo> {
	static const EMap<MaskingAlgo> to_string;
	static const SEMap<MaskingAlgo> from_string;
};

enum class MaskingMode { NONE, TANTAN, BLAST_SEG };

template<>
struct EnumTraits<MaskingMode> {
	static const SEMap<MaskingMode> from_string;
};

struct MaskingTable {

	MaskingTable();
	MaskingTable(const MaskingTable& t);
	MaskingTable& operator=(const MaskingTable& t);
	void add(const size_t block_id, const int begin, const int end, Letter* seq);
	void remove(SequenceSet& seqs, const int template_len, const bool add_bit_mask) const;
	void apply(SequenceSet& seqs) const;
	bool blank() const;
	size_t masked_letters() const;

private:

	struct Entry {
		Entry(const size_t block_id, const int begin) :
			block_id(block_id), begin(begin) {}
		size_t block_id;
		int begin;
	};

	size_t seq_count_, masked_letters_;
	std::vector<Entry> entry_;
	StringSetBase<Letter, Sequence::DELIMITER, 1> seqs_;
	std::mutex mtx_;

};

template<size_t K>
struct Kmer {
	Kmer():
		code(0)
	{ }
	Kmer(const char* s) :
		code(0)
	{
		assert(strlen(s) == K);
		for (size_t i = 0; i < K; ++i)
			code = (code << 5) | (uint64_t)amino_acid_traits.from_char(*s++);
	}
	operator uint64_t() const {
		return code;
	}
	uint64_t code;
};

template<size_t K>
struct std::hash<Kmer<K>> {
	size_t operator()(const Kmer<K> k) const {
		return std::hash<uint64_t>()(k.code);
	}
};

template<size_t K>
struct KmerIterator {
	KmerIterator(const Sequence& seq) :
		ptr_(seq.data() - 1),
		end_(seq.end())
	{
		inc(0);
	}
	Kmer<K> operator*() const {
		return kmer_;
	}
	bool good() const {
		return ptr_ < end_;
	}
	KmerIterator& operator++() {
		inc(K - 1);
		return *this;
	}
	ptrdiff_t operator-(const Letter* ptr) {
		return ptr_ + 1 - K - ptr;
	}
private:
	static const uint64_t BITS = 5;
	void inc(size_t n) {
		do {
			++ptr_;
			if (ptr_ == end_)
				return;
			const uint64_t l = letter_mask(*ptr_);
			if (l < TRUE_AA) {
				kmer_.code = (kmer_.code << BITS) | l;
				++n;
			}
			else
				n = 0;
		} while (n < K);
		kmer_.code &= (uint64_t(1) << (BITS * K)) - 1;
	}
	const Letter* ptr_, *const end_;
	Kmer<K> kmer_;
};

const size_t MOTIF_LEN = 8;
extern const std::unordered_set<Kmer<MOTIF_LEN>> motif_table;