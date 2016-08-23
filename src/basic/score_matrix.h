/****
Copyright (c) 2016, University of Tuebingen, Benjamin Buchfink
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

#ifndef SCORE_MATRIX_H_
#define SCORE_MATRIX_H_

#include <limits>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include "../util/log_stream.h"
#include "value.h"

using std::string;
using std::cout;
using std::endl;

const double LN_2 = 0.69314718055994530941723212145818;

struct Score_matrix
{

	Score_matrix() {}
	Score_matrix(const string &matrix, int gap_open, int gap_extend, int reward, int penalty);
	Score_matrix(const string &matrix_file, double lambda, double K, int gap_open, int gap_extend);

	friend std::ostream& operator<<(std::ostream& s, const Score_matrix &m)
	{
		s << "(Matrix=" << m.name_
			<< " Lambda=" << m.lambda()
			<< " K=" << m.k()
			<< " Penalties=" << m.gap_open_
			<< '/' << m.gap_extend_ << ')';
		return s;
	}

	const int8_t* matrix8() const
	{ return matrix8_.data; }

	const uint8_t* matrix8u() const
	{ return matrix8u_.data; }

	const int16_t* matrix16() const
	{ return matrix16_.data; }

	int operator()(Letter a, Letter b) const
	{
		return matrix8_.data[(int(a) << 5) + int(b)];
	}

	uint8_t biased_score(Letter a, Letter b) const
	{ return matrix8u_.data[(int(a) << 5) + int(b)]; }

	char bias() const
	{ return bias_; }

	double bitscore(int raw_score) const
	{ return ( lambda() * raw_score - ln_k()) / LN_2; }

	double rawscore(double bitscore, double) const
	{ return (bitscore*LN_2 + ln_k()) / lambda(); }

	int rawscore(double bitscore) const
	{ return (int)ceil(rawscore(bitscore, double ())); }

	double evalue(int raw_score, size_t db_letters, unsigned query_len) const
	{ return static_cast<double>(db_letters) * query_len * pow(2,-bitscore(raw_score)); }

	double bitscore(double evalue, size_t db_letters, unsigned query_len) const
	{ return -log(evalue/db_letters/query_len)/log(2); }

	double lambda() const
	{
		return constants_[3];
	}

	double k() const
	{
		return constants_[4];
	}

	double ln_k() const
	{
		return log(k());
	}

	char low_score() const;

private:

	template<typename _t>
	struct Scores
	{
		Scores() {}
		Scores(const char *scores, char bias = 0)
		{
			const unsigned n = value_traits.alphabet_size;
			for(unsigned i=0;i<32;++i)
				for(unsigned j=0;j<32;++j)
					data[i*32+j] = i < n && j < n ? (_t)(scores[i*n+j] + (int)bias) : std::numeric_limits<_t>::min();
		}
#ifdef _MSC_VER
		__declspec(align(16)) _t data[32 * 32];
#else
		_t data[32 * 32] __attribute__((aligned(16)));
#endif
	};

	int gap_open_, gap_extend_;
	const double *constants_;
	string name_;
	Scores<int8_t> matrix8_;
	char bias_;
	Scores<uint8_t> matrix8u_;
	Scores<int16_t> matrix16_;

};

extern Score_matrix score_matrix;

#endif /* SCORE_MATRIX_H_ */
