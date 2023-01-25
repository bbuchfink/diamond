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

#pragma once
#include <limits>
#include <ostream>
#include <math.h>
#include <stdint.h>
#include "../util/log_stream.h"
#include "../basic/value.h"
#include "../lib/alp/sls_alignment_evaluer.hpp"
#include "../stats/standard_matrix.h"

const double LN_2 = 0.69314718055994530941723212145818;

struct ScoreMatrix
{

	struct Custom {};

	ScoreMatrix() :ln_k_(0.0) {}
	ScoreMatrix(const std::string& matrix, int gap_open, int gap_extend, int frame_shift, int stop_match_score, uint64_t db_letters = 0, int scale = 1, bool mmseqs_compat = false);
	ScoreMatrix(const std::string &matrix_file, int gap_open, int gap_extend, int stop_match_score, const Custom&, uint64_t db_letters = 0);

	friend std::ostream& operator<<(std::ostream& s, const ScoreMatrix&m);

	const int8_t* matrix8() const
	{ return matrix8_.data; }

	const int8_t* matrix8_low() const
	{
		return matrix8_low_.data;
	}

	const int8_t* matrix8_high() const
	{
		return matrix8_high_.data;
	}

	const int8_t* matrix8u_low() const
	{
		return matrix8u_low_.data;
	}

	const int8_t* matrix8u_high() const
	{
		return matrix8u_high_.data;
	}

	const uint8_t* matrix8u() const
	{ return matrix8u_.data; }

	const int16_t* matrix16() const
	{ return matrix16_.data; }

	const int* matrix32() const
	{
		return matrix32_.data;
	}

	std::vector<const int*> matrix32_scaled_pointers() const {
		return matrix32_scaled_.pointers();
	}

	int operator()(Letter a, Letter b) const
	{
		return matrix32_.data[(int(a) << 5) + int(b)];
	}

	int operator()(size_t a, size_t b) const
	{
		return matrix32_.data[(a << 5) + b];
	}

	const int* row(Letter a) const
	{
		return &matrix32_.data[(int)a << 5];
	}

	uint8_t biased_score(Letter a, Letter b) const
	{ return matrix8u_.data[(int(a) << 5) + int(b)]; }

	char bias() const
	{ return bias_; }

	double bitscore(double raw_score) const;

	double rawscore(double bitscore, double) const
	{ return (bitscore*LN_2 + ln_k()) / lambda(); }

	int rawscore(double bitscore) const
	{ return (int)ceil(rawscore(bitscore, double ())); }

	double evalue(int raw_score, unsigned query_len, unsigned subject_len) const;
	double evalue_norm(int raw_score, unsigned query_len, unsigned subject_len) const;
	double bitscore_corrected(int raw_score, unsigned query_len, unsigned subject_len) const;

	double evalue_norm(int raw_score, int query_len) const
	{
		return 1e9 * query_len * pow(2, -bitscore(raw_score * scale_));
	}

	/*double bitscore(double evalue, unsigned query_len) const
	{ return -log(evalue/db_letters_/query_len)/log(2); }*/

	double bitscore_norm(double evalue, unsigned query_len) const
	{
		return -log(evalue / 1e9 / query_len) / log(2);
	}

	double lambda() const
	{
		return evaluer.parameters().lambda;
	}

	double k() const
	{
		return evaluer.parameters().K;
	}

	double ln_k() const
	{
		return ln_k_;
	}

	int8_t low_score() const;
	int8_t high_score() const;

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

	const double* joint_probs() const {
		return (const double*)standard_matrix_->joint_probs;
	}

	const double* background_freqs() const {
		return standard_matrix_->background_freqs.data();
	}

	double ungapped_lambda() const {
		return standard_matrix_->ungapped_constants().Lambda;
	}

	double ideal_lambda() const {
		return ideal_lambda_;
	}

	const Stats::FreqRatios& freq_ratios() const {
		return standard_matrix_->freq_ratios;
	}

	const std::array<double, TRUE_AA>& background_scores() const {
		return background_scores_;
	}

	std::string name() const {
		return name_;
	}

	double avg_id_score() const;
	bool report_cutoff(int score, double evalue) const;

private:

	template<typename _t>
	struct Scores
	{
		Scores() {}
		Scores(const int8_t *scores, int stop_match_score = 1, int8_t bias = 0, unsigned modulo = 32, unsigned offset = 0)
		{
			const unsigned n = value_traits.alphabet_size;
			for (unsigned i = 0; i < 32; ++i)
				for (unsigned j = 0; j < 32; ++j) {
					const unsigned j2 = j % modulo + offset;
					data[i * 32 + j] = i < n && j2 < n ? (_t)(scores[i * n + j2] + (int)bias) : SCHAR_MIN;
				}
			if (stop_match_score != 1)
				data[24 * 32 + 24] = stop_match_score;
		}
		Scores(const double (&freq_ratios)[Stats::NCBI_ALPH][Stats::NCBI_ALPH], double lambda, const int8_t* scores, int scale);
		std::vector<const _t*> pointers() const;
		friend std::ostream& operator<<(std::ostream& s, const Scores& scores) {
			for (int i = 0; i < 20; ++i) {
				for (int j = 0; j < 20; ++j)
					s << scores.data[i * 32 + j] << '\t';
				s << std::endl;
			}
			return s;
		}
		alignas(32) _t data[32 * 32];
	};

	void init_background_scores();

	const Stats::StandardMatrix* standard_matrix_;
	const int8_t* score_array_;
	int gap_open_, gap_extend_, frame_shift_;
	double db_letters_, ln_k_, scale_;
	std::string name_;
	Scores<int8_t> matrix8_;
	Scores<int> matrix32_, matrix32_scaled_;
	int8_t bias_;
	double ideal_lambda_;
	Scores<uint8_t> matrix8u_;
	Scores<int8_t> matrix8_low_;
	Scores<int8_t> matrix8_high_;
	Scores<int8_t> matrix8u_low_;
	Scores<int8_t> matrix8u_high_;
	Scores<int16_t> matrix16_;
	std::array<double, TRUE_AA> background_scores_;
	Sls::AlignmentEvaluer evaluer;

};

extern ScoreMatrix score_matrix;
