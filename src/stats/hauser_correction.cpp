/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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

#include <array>
#include <algorithm>
#include "../dp/dp.h"
#include "score_matrix.h"

using std::array;
using std::vector;

static const int PADDING = 32;

struct Vector_scores
{
	Vector_scores()
	{
		memset(scores, 0, sizeof(scores));
	}
	Vector_scores& operator+=(Letter l)
	{
		for (unsigned i = 0; i < 20; ++i)
			scores[i] += score_matrix(l, i);
		return *this;
	}
	Vector_scores& operator-=(Letter l)
	{
		for (unsigned i = 0; i < 20; ++i)
			scores[i] -= score_matrix(l, i);
		return *this;
	}
	int scores[20];
};

Bias_correction::Bias_correction(const Sequence &seq):
	vector<float>(seq.length())
{
	Vector_scores scores;
	const unsigned window = config.cbs_window, window_half = std::min(window/2, (unsigned)seq.length() - 1);
	unsigned n = 0;
	unsigned h = 0, m = 0, t = 0, l = (unsigned)seq.length();
	const std::array<double, TRUE_AA>& background_scores = score_matrix.background_scores();
	while (n < window_half && h < l) {
		++n;
		scores += seq[h];
		++h;
	}
	while (n < (window+1) && h < l) {
		++n;
		scores += seq[h];
		const Letter r = seq[m];
		if (r < 20)
			//this->operator[](m) = std::min((float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1), 0.0f);
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++h;
		++m;
	}
	while (h < l) {
		scores += seq[h];
		scores -= seq[t];
		const Letter r = seq[m];
		if (r < 20)
			//this->operator[](m) = std::min((float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1), 0.0f);
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++h;
		++t;
		++m;
	}
	while (m < l && n>(window_half+1)) {
		--n;
		scores -= seq[t];
		const Letter r = seq[m];
		if (r < 20)
			//this->operator[](m) = std::min((float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1), 0.0f);
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++t;
		++m;
	}
	while (m < l) {
		const Letter r = seq[m];
		if (r < 20)
			//this->operator[](m) = std::min((float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1), 0.0f);
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++m;
	}

	int8.reserve(seq.length() + PADDING);
	for (float f : *this)
		int8.push_back(int8_t(f < 0.0f ? f - 0.5f : f + 0.5f));
	int8.insert(int8.end(), PADDING, 0);
}

int Bias_correction::operator()(const Hsp &hsp) const
{
	float s = 0;
	for (Hsp::Iterator i = hsp.begin(); i.good(); ++i) {
		switch (i.op()) {
		case op_match:
		case op_substitution:
			s += (*this)[i.query_pos.translated];
		default:
			;
		}
	}
	return (int)s;
}

int Bias_correction::operator()(const DiagonalSegment &d) const
{
	float s = 0;
	const int end = d.query_end();
	for (int i = d.i; i < end; ++i)
		s += (*this)[i];
	return (int)s;
}

vector<int8_t> Bias_correction::reverse(const int8_t* p, const size_t len) {
	vector<int8_t> r;
	if (p == nullptr)
		return r;
	r.reserve(len);
	std::reverse_copy(p, p + len, std::back_inserter(r));
	return r;
}

namespace Stats {

vector<int> hauser_global(const Composition& query_comp, const Composition& target_comp) {
	const std::array<double, TRUE_AA>& background_scores = score_matrix.background_scores();
	array<double, TRUE_AA> qscores, tscores;
	qscores.fill(0.0);
	tscores.fill(0.0);
	for (size_t i = 0; i < TRUE_AA; ++i)
		for (size_t j = 0; j < TRUE_AA; ++j) {
			qscores[i] += query_comp[j] * (double)score_matrix(i, j);
			tscores[i] += target_comp[j] * (double)score_matrix(i, j);
		}

	for (size_t i = 0; i < TRUE_AA; ++i) {
		qscores[i] = (background_scores[i] - qscores[i]);
		tscores[i] = (background_scores[i] - tscores[i]);
	}

	vector<int> m(AMINO_ACID_COUNT * AMINO_ACID_COUNT);
	for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
		for (size_t j = 0; j < AMINO_ACID_COUNT; ++j) {
			double s = (double)score_matrix(i, j), q = i < TRUE_AA ? qscores[i] : 0.0, t = j < TRUE_AA ? tscores[j] : 0.0;
			m[i * AMINO_ACID_COUNT + j] = (int)std::round(s + std::min(q, t));
		}
	return m;
}

}