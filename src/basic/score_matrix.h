/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef SCORE_MATRIX_H_
#define SCORE_MATRIX_H_

#include "../algo/blast/core/blast_stat.h"
#include "../algo/blast/core/blast_encoding.h"
#include "score_traits.h"

using std::string;
using std::cout;
using std::endl;
using std::auto_ptr;

struct Score_params_exception : public std::exception
{
	virtual const char* what() const throw()
	{ return "Invalid scoring parameters"; }
};

struct Blast_score_blk
{

	template<typename _val>
	Blast_score_blk(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const _val&):
		data_ (BlastScoreBlkNew(blast_seq_code<_val>(), 1))
	{
		if(data_ == 0)
			throw Score_params_exception ();
		if((data_->kbp_gap_std[0] = Blast_KarlinBlkNew()) == 0
				|| (data_->kbp_std[0] = Blast_KarlinBlkNew()) == 0)
			throw Score_params_exception ();
		if(blast_load_karlin_blk<_val>(data_->kbp_gap_std[0],
				data_->kbp_std[0],
				gap_open,
				gap_extend,
				reward,
				penalty,
				matrix.c_str()) != 0)
			throw Score_params_exception ();
		data_->name = const_cast<char*>(matrix.c_str());
		data_->reward = reward;
		data_->penalty = penalty;
		if(Blast_ScoreBlkMatrixFill (data_, 0) != 0)
			throw Score_params_exception ();
		data_->name = 0;
	}

	~Blast_score_blk()
	{ BlastScoreBlkFree(data_); }

	template<typename _val>
	int score(_val x, _val y) const
	{ return data_->matrix->data[(long)blast_alphabet<_val>()[(long)Value_traits<_val>::ALPHABET[x]]][(long)blast_alphabet<_val>()[(long)Value_traits<_val>::ALPHABET[y]]]; }

	double lambda() const
	{ return data_->kbp_gap_std[0]->Lambda; }

	double k() const
	{ return data_->kbp_gap_std[0]->K; }

	double ln_k() const
	{ return data_->kbp_gap_std[0]->logK; }

	int low_score() const
	{ return data_->loscore; }

private:

	BlastScoreBlk *data_;

};

struct score_matrix
{

	template<typename _val>
	score_matrix(const string &matrix, int gap_open, int gap_extend, int reward, int penalty, const _val&):
		sb_ (matrix, gap_open, gap_extend, reward, penalty, _val ()),
		bias_ ((char)(-sb_.low_score())),
		name_ (matrix),
		matrix8_ (_val(), sb_),
		matrix8u_ (_val(), sb_, bias_),
		matrix16_ (_val(), sb_)
	{ }

	template<typename _val>
	void print() const
	{
		verbose_stream << "Scoring matrix = " << name_ << endl;
		verbose_stream << "Lambda = " << sb_.lambda() << endl;
		verbose_stream << "K = " << sb_.k() << endl;
		/*const unsigned n = Value_traits<_val>::ALPHABET_SIZE;
		for(unsigned i=0;i<n;++i) {
			for(unsigned j=0;j<n;++j)
				printf("%3i", (int)matrix8_.data[i*32+j]);
			printf("\n");
		}*/
	}

	static const score_matrix& get()
	{ return *instance; }

	const int8_t* matrix8() const
	{ return matrix8_.data; }

	const uint8_t* matrix8u() const
	{ return matrix8u_.data; }

	const int16_t* matrix16() const
	{ return matrix16_.data; }

	template<typename _val>
	int letter_score(_val a, _val b) const
	{ return matrix8_.data[(int(a) << 5) + int(b)]; }

	template<typename _val>
	uint8_t biased_score(_val a, _val b) const
	{ return matrix8u_.data[(int(a) << 5) + int(b)]; }

	char bias() const
	{ return bias_; }

	double bitscore(int raw_score) const
	{ return ( sb_.lambda() * raw_score - sb_.ln_k()) / LN_2; }

	double rawscore(double bitscore, double) const
	{ return (bitscore*LN_2 + sb_.ln_k()) / sb_.lambda(); }

	int rawscore(double bitscore) const
	{ return (int)ceil(rawscore(bitscore, double ())); }

	double evalue(int raw_score, size_t db_letters, unsigned query_len) const
	{ return static_cast<double>(db_letters) * query_len * pow(2,-bitscore(raw_score)); }

	double bitscore(double evalue, size_t db_letters, unsigned query_len) const
	{ return -log(evalue/db_letters/query_len)/log(2); }

	double k() const
	{ return sb_.k(); }

	double lambda() const
	{ return sb_.lambda(); }

	static auto_ptr<score_matrix> instance;

private:

	static const double	LN_2 = 0.69314718055994530941723212145818;

	template<typename _t>
	struct Scores
	{
		template<typename _val>
		Scores(const _val&, const Blast_score_blk &sb, char bias = 0)
		{
			const unsigned n = Value_traits<_val>::ALPHABET_SIZE;
			for(unsigned i=0;i<32;++i)
				for(unsigned j=0;j<32;++j)
					data[i*32+j] = i < n && j < n ? (_t)(sb.score((_val)i, (_val)j) + (int)bias) : std::numeric_limits<_t>::min();
		}
		_t data[32*32] __attribute__ ((aligned (16)));
	};

	const Blast_score_blk sb_;
	const char bias_;
	const string name_;
	const Scores<int8_t> matrix8_;
	const Scores<uint8_t> matrix8u_;
	const Scores<int16_t> matrix16_;

};

auto_ptr<score_matrix> score_matrix::instance;

#endif /* SCORE_MATRIX_H_ */
