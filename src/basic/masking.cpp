/****
Copyright (c) 2017, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include "masking.h"

Masking::Masking(const Score_matrix &score_matrix)
{
	const double lambda = 0.324032;//  score_matrix.lambda();		//;
	for (int i = 0; i < size; ++i) {
		mask_table_[i] = value_traits.mask_char;
		for (int j = 0; j < size; ++j)
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
		0.5, (const tantan::uchar*)mask_table_);
}

void mask_worker(Atomic<size_t> *next, Sequence_set *seqs, const Masking *masking)
{
	size_t i;
	while ((i = (*next)++) < seqs->get_length())
		masking->operator()(seqs->ptr(i), seqs->length(i));
}

void mask_seqs(Sequence_set &seqs, const Masking &masking)
{
	Thread_pool threads;
	Atomic<size_t> next(0);
	for (size_t i = 0; i < config.threads_; ++i)
		threads.push_back(launch_thread(mask_worker, &next, &seqs, &masking));
	threads.join_all();
}
