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

#include <math.h>
#include "masking.h"
#include "../lib/tantan/tantan.hh"

using namespace std;

unique_ptr<Masking> Masking::instance;
const uint8_t Masking::bit_mask = 128;

Masking::Masking(const Score_matrix &score_matrix)
{
	const double lambda = score_matrix.lambda(); // 0.324032
	for (unsigned i = 0; i < size; ++i) {
		mask_table_x_[i] = value_traits.mask_char;
		mask_table_bit_[i] = (uint8_t)i | bit_mask;
		for (unsigned j = 0; j < size; ++j)
			if (i < value_traits.alphabet_size && j < value_traits.alphabet_size)
				likelihoodRatioMatrix_[i][j] = exp(lambda * score_matrix(i, j));
	}
	std::copy(likelihoodRatioMatrix_, likelihoodRatioMatrix_ + size, probMatrixPointers_);
	int firstGapCost = score_matrix.gap_extend() + score_matrix.gap_open();
	firstGapProb_ = exp(-lambda * firstGapCost);
	otherGapProb_ = exp(-lambda * score_matrix.gap_extend());
	firstGapProb_ /= (1 - otherGapProb_);
}

void Masking::operator()(Letter *seq, size_t len) const
{
	tantan::maskSequences((tantan::uchar*)seq, (tantan::uchar*)(seq + len), 50,
		(tantan::const_double_ptr*)probMatrixPointers_,
		0.005, 0.05,
		0.9,
		0, 0,
		0.5, (const tantan::uchar*)mask_table_x_);
}

void Masking::mask_bit(Letter *seq, size_t len) const
{
	tantan::maskSequences((tantan::uchar*)seq, (tantan::uchar*)(seq + len), 50,
		(tantan::const_double_ptr*)probMatrixPointers_,
		0.005, 0.05,
		0.9,
		0, 0,
		0.5, (const tantan::uchar*)mask_table_bit_);
}

void Masking::bit_to_hard_mask(Letter *seq, size_t len, size_t &n) const
{
	for (size_t i = 0; i < len; ++i)
		if (seq[i] & bit_mask) {
			seq[i] = value_traits.mask_char;
			++n;
		}
}

void Masking::remove_bit_mask(Letter *seq, size_t len) const
{
	for (size_t i = 0; i < len; ++i)
		if (seq[i] & bit_mask)
			seq[i] &= ~bit_mask;
}

void mask_worker(Atomic<size_t> *next, Sequence_set *seqs, const Masking *masking, bool hard_mask)
{
	size_t i;
	while ((i = (*next)++) < seqs->get_length())
		if (hard_mask)
			masking->operator()(seqs->ptr(i), seqs->length(i));
		else
			masking->mask_bit(seqs->ptr(i), seqs->length(i));
}

void mask_seqs(Sequence_set &seqs, const Masking &masking, bool hard_mask)
{
	Thread_pool threads;
	Atomic<size_t> next(0);
	for (size_t i = 0; i < config.threads_; ++i)
		threads.push_back(launch_thread(mask_worker, &next, &seqs, &masking, hard_mask));
	threads.join_all();
}