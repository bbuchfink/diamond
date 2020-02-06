/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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
	Score_matrix(const string &matrix, int gap_open, int gap_extend, int frame_shift, int stop_match_score, uint64_t db_letters = 0);
	Score_matrix(const string &matrix_file, double lambda, double K, int gap_open, int gap_extend, uint64_t db_letters = 0);

	friend std::ostream& operator<<(std::ostream& s, const Score_matrix &m);

	const int8_t* matrix8() const
	{ return matrix8_.data; }

	const uint8_t* matrix8u() const
	{ return matrix8u_.data; }

	const int16_t* matrix16() const
	{ return matrix16_.data; }

	const int* matrix32() const
	{
		return matrix32_.data;
	}

	int operator()(Letter a, Letter b) const
	{
		return matrix8_.data[(int(a) << 5) + int(b)];
	}

	const int* row(char a) const
	{
		return &matrix32_.data[(int)a << 5];
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

	double evalue(int raw_score, unsigned query_len) const
	{ return db_letters_ * query_len * pow(2, -bitscore(raw_score)); }

	double bitscore(double evalue, unsigned query_len) const
	{ return -log(evalue/db_letters_/query_len)/log(2); }

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
	char high_score() const;

	int gap_open() const
	{
		return gap_open_;
	}

	int gap_extend() const
	{
		return gap_extend_;
	}

	int frame_shift() const
	{
		return frame_shift_;
	}

	uint64_t db_letters() const
	{
		return (uint64_t)db_letters_;
	}

	void set_db_letters(uint64_t n)
	{
		db_letters_ = (double)n;
	}

	double avg_id_score() const;

private:

	template<typename _t>
	struct Scores
	{
		Scores() {}
		Scores(const char *scores, int stop_match_score = 1, char bias = 0)
		{
			const unsigned n = value_traits.alphabet_size;
			for (unsigned i = 0; i < 32; ++i)
				for (unsigned j = 0; j < 32; ++j)
					data[i * 32 + j] = i < n && j < n ? (_t)(scores[i*n + j] + (int)bias) : -(std::numeric_limits<_t>::max() / 2);
			if (stop_match_score != 1)
				data[24 * 32 + 24] = stop_match_score;
		}
#ifdef _MSC_VER
		__declspec(align(16)) _t data[32 * 32];
#else
		_t data[32 * 32] __attribute__((aligned(16)));
#endif
	};

	int gap_open_, gap_extend_, frame_shift_;
	double db_letters_;
	const double *constants_;
	string name_;
	Scores<int8_t> matrix8_;
	char bias_;
	Scores<uint8_t> matrix8u_;
	Scores<int16_t> matrix16_;
	Scores<int> matrix32_;

};

extern Score_matrix score_matrix;
typedef signed char MatrixTable[AMINO_ACID_COUNT*AMINO_ACID_COUNT];
extern const MatrixTable s_Blosum45PSM, s_Blosum50PSM, s_Blosum62PSM, s_Blosum80PSM, s_Blosum90PSM, s_Pam250PSM, s_Pam30PSM, s_Pam70PSM;

#endif /* SCORE_MATRIX_H_ */
