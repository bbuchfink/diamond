/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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
#include "value.h"
#include "../stats/score_matrix.h"
#include "../basic/sequence.h"
#include "../data/sequence_set.h"
#include "../lib/blast/blast_seg.h"

struct Masking
{
	enum struct Algo { TANTAN, SEG };

	Masking(const Score_matrix &score_matrix);
	~Masking();
	void operator()(Letter *seq, size_t len, Algo algo = Algo::TANTAN) const;
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

size_t mask_seqs(SequenceSet &seqs, const Masking &masking, bool hard_mask = true, Masking::Algo algo = Masking::Algo::TANTAN);
