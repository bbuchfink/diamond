/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include <vector>
#include <math.h>
#include "value.h"
#include "score_matrix.h"
#include "../basic/sequence.h"
#include "../lib/tantan/tantan.hh"
#include "../data/sequence_set.h"

using std::vector;

struct Masking
{
	Masking(const Score_matrix &score_matrix);
	void operator()(Letter *seq, size_t len) const;
	void mask_bit(Letter *seq, size_t len) const;
	void bit_to_hard_mask(Letter *seq, size_t len, size_t &n) const;
	void remove_bit_mask(Letter *seq, size_t len) const;
	static const Masking& get()
	{
		return *instance;
	}
	static auto_ptr<Masking> instance;
	static const uint8_t bit_mask;
private:
	enum { size = 64 };
	double likelihoodRatioMatrix_[size][size], *probMatrixPointers_[size], firstGapProb_, otherGapProb_;
	char mask_table_x_[size], mask_table_bit_[size];
};

void mask_seqs(Sequence_set &seqs, const Masking &masking, bool hard_mask = true);