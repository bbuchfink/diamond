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

#include "dp.h"

int Diag_scores::min_diag_score = 19;
int Diag_scores::min_low_score = 13;

void score_diagonal2(const Letter *query, const Bias_correction &query_bc, const Letter *subject, int len, int qbegin, int jbegin, vector<Diagonal_node> &diags, int cutoff)
{
	typedef int score_t;
	static const score_t xdrop = 10;
	int i = 0, j = 0, begin = 0, end = 0, l = 0;
	score_t score = 0, max_score = 0;
	while (l<len) {
		score += score_matrix(query[i], subject[i]); // +query_bc[qbegin + i];
		if (score <= 0 || max_score - score > xdrop) {
			if (max_score >= cutoff) {
				diags.push_back(Diagonal_node(qbegin + begin, jbegin + begin, end - begin, (int)max_score));
				//cells += end - begin;
			}
			score = 0;
			max_score = 0;
			begin = i + 1;
		}
		else if (score > max_score) {
			max_score = score;
			end = i + 1;
		}
		++i;
		++l;
	}
	if (max_score >= cutoff) {
		diags.push_back(Diagonal_node(qbegin + begin, jbegin + begin, end - begin, (int)max_score));
		//cells += end - begin;
	}
}

void scan_cols(const Long_score_profile &qp, sequence s, int i, int j, int j_end, vector<uint8_t> &sv_max, bool log, Band &buf, Band &local_max, int block_len)
{
#ifdef __SSE2__
	typedef score_vector<uint8_t> Sv;
	const Sv vbias(score_matrix.bias());
	const int qlen = (int)qp.length(),
		diags = buf.diags();

	int j2 = std::max(-(i - j + 15), j),
		i3 = j2 + i - j,
		j2_end = std::min(qlen - (i - j), j_end);
	uint8_t *local_max_ptr = local_max.data() + (j2 - j) / 16 * diags,
		*buf_ptr = buf.data() + (j2 - j)*diags;
	Sv v, max, global_max;
	for (; j2 < j2_end; ++j2, ++i3) {
		assert(j2 >= 0);
		const uint8_t *q = qp.get(s[j2], i3);
		v = v + score_vector<uint8_t>(q);
		v -= vbias;
		max.max(v);
		assert(buf.check(buf_ptr + 16));
		v.store(buf_ptr);
		buf_ptr += diags;
		if (((j2 - j) & 15) == 15) {
			global_max.max(max);
			assert(local_max.check(local_max_ptr + 16));
			max.store(local_max_ptr);
			local_max_ptr += diags;
			max = score_vector<uint8_t>();
		}
	}
	if (((j2 - j) & 15) != 0) {
		global_max.max(max);
		assert(local_max.check(local_max_ptr));
		max.store(local_max_ptr);
	}
	global_max.store(&sv_max[0]);
#endif
}

void scan_cols(const sequence &q, sequence s, int i, int j, int j_end, vector<uint8_t> &sv_max, bool log, Band &buf, Band &local_max, int block_len)
{
#ifdef __SSE2__
	typedef score_vector<uint8_t> Sv;
	const Sv vbias(score_matrix.bias());
	const int qlen = (int)q.length(),
		diags = buf.diags();

	int j2 = std::max(-(i - j + 15), j),
		i3 = j2 + i - j,
		j2_end = std::min(qlen - (i - j), j_end);
	uint8_t *local_max_ptr = local_max.data() + (j2 - j) / 16 * diags,
		*buf_ptr = buf.data() + (j2 - j)*diags;
	Sv v, max, global_max;
	//__m128i f = _mm_loadu_si128((__m128i*)&q[0]);
	for (; j2 < j2_end; ++j2, ++i3) {
		assert(j2 >= 0);
		const Sv scores(s[j2], _mm_loadu_si128((__m128i*)&q[i3]));
		//const Sv scores(s[j2], f);
		v = v + scores;
		v -= vbias;
		max.max(v);
		assert(buf.check(buf_ptr + 16));
		v.store(buf_ptr);
		buf_ptr += diags;
		if (((j2 - j) & 15) == 15) {
			global_max.max(max);
			assert(local_max.check(local_max_ptr + 16));
			max.store(local_max_ptr);
			local_max_ptr += diags;
			max = score_vector<uint8_t>();
		}
	}
	if (((j2 - j) & 15) != 0) {
		global_max.max(max);
		assert(local_max.check(local_max_ptr));
		max.store(local_max_ptr);
	}
	global_max.store(&sv_max[0]);
#else
	const int qlen = (int)q.length(),
		diags = buf.diags();

	for (int i0 = i; i0 < i + 16; ++i0) {
		int j2 = std::max(-(i0 - j), j),
			i3 = j2 + i0 - j,
			j2_end = std::min(qlen - (i0 - j), j_end);
		uint8_t *local_max_ptr = local_max.data() + (j2 - j) / 16 * diags + (i0 - i),
			*buf_ptr = buf.data() + (j2 - j)*diags + (i0 - i);
		int v = 0, max = 0, global_max = 0;
		for (; j2 < j2_end; ++j2, ++i3) {
			assert(j2 >= 0);
			v = std::max(v + score_matrix(q[i3], s[j2]), 0);
			max = std::max(max, v);
			assert(buf.check(buf_ptr + 1));
			*buf_ptr = (uint8_t)std::min(v, 255);
			buf_ptr += diags;
			if (((j2 - j) & 15) == 15) {
				global_max = std::max(global_max, max);
				assert(local_max.check(local_max_ptr + 1));
				*local_max_ptr = (uint8_t)std::min(max,255);
				local_max_ptr += diags;
				max = 0;
			}
		}
		if (((j2 - j) & 15) != 0) {
			global_max = std::max(global_max, max);
			assert(local_max.check(local_max_ptr));
			*local_max_ptr = (uint8_t)std::min(max, 255);
		}
		sv_max[i0 - i] = (uint8_t)std::min(global_max, 255);
	}
#endif
}


int get_low_idx(Band::Iterator &d, int begin, int end, int d0)
{
	uint8_t low = 255;
	int j = end;
	for (int i = end - 1; i >= begin; --i)
		if (d[i] == 0)
			return i;
		else if (d[i] < low) {
			low = d[i];
			j = i;
		}
	return begin > d0 ? j : d0 - 1;
}

int get_max_idx(Band::Iterator &d, int begin, int end)
{
	assert(begin >= 0 && begin < end);
	int i = begin, s = d[begin];
	for (int j = i + 1; j < end; ++j)
		if (d[j] > s) {
			i = j;
			s = d[j];
		}
	return i;
}

int get_diag(int i, int j, Band::Iterator &d, int begin, int last, int end, int d0, vector<Diagonal_node> &diags, int block_len, bool log, int cutoff, int best_score, const Bias_correction &query_bc)
{
	assert(end >= begin && begin >= 0);
	int z = std::numeric_limits<int>::max();
	while (end > begin) {
		/*const int p1 = get_score_idx(d, std::max(end - block_len, begin), end, max_score),
		p0 = get_low_idx(d, last, std::min(std::min(begin + block_len, end), p1));*/
		const int mod = end%block_len,
			p1 = get_max_idx(d, std::max(begin, end - (mod == 0 ? block_len : mod)), end),
			p0 = get_low_idx(d, last, p1, d0);
		assert(p1 >= p0);
		assert(p1 < end);
		const int score = d[p1] - (p0 >= d0 ? d[p0] : 0);
		if (score >= cutoff || (d[p1] == best_score &&  score > 0)) {
			//if (score >= cutoff) {
			assert(i + p0 + 1 >= 0);
			assert(j + p0 + 1 >= 0);
			Diagonal_segment diag(Diagonal_segment(i + p0 + 1, j + p0 + 1, p1 - p0, score));
			//if (diag.score + query_bc(diag) >= cutoff) {
			diag.score = diag.score; // +query_bc(diag);
			diags.push_back(diag);
			assert(p0 + 1 >= 0);
			z = p0 + 1;
			/*if (log)
			cout << diags.back() << endl;*/
			//}
		}
		end = p0;
	}
	assert(z >= last);
	return z;
}

void Diag_scores::get_diag(int i, int j, int o, int j_begin, int j_end, vector<Diagonal_node> &diags, int cutoff, bool log, const Bias_correction &query_bc)
{
	Band::Iterator d(local_max.diag(o)), d2(score_buf.diag(o));
	const int diag = i - j,
		j0 = dj0(diag),
		j1 = dj1(diag),
		b0 = (j0 - j) / block_len,
		b1 = (j1 - j + block_len - 1) / block_len;
	int p = std::max((j_begin - j) / block_len, b0), p_begin = p, p_end = std::min((j_end - j + block_len - 1) / block_len, b1), best = -1, best_score = -1, begin = -1, max_score, last = p;
	while (last > b0)
		--last;
	for (; p < p_end; ++p) {
		if ((max_score = d[p]) >= cutoff
			&& (p == 0 || max_score > d[p - 1])) {
			if (begin == -1)
				begin = p;
			best = p;
			best_score = max_score;
		}
		else if (begin != -1) {
			const int z = ::get_diag(i, j, d2, std::max(begin*block_len, j0 - j), std::max(last*block_len, j0 - j), std::min((best + 1)*block_len, j1 - j), j0 - j, diags, block_len, log, cutoff, best_score, query_bc);
			if (z < std::numeric_limits<int>::max()) {
				assert(diags.back().len > 0);
				assert(diags.back().j >= 0 && diags.back().subject_end() <= slen);
				assert(diags.back().i >= 0 && diags.back().query_end() <= qlen);
				last = best + 1;
			}
			begin = -1;
			best = -1;
		}
	}

	if (begin != -1) {
		if (best == p_end - 1) {
			for (; best<b1 && d[best] >= cutoff && (best == 0 || d[best] > d[best - 1]); ++best);
			best -= 1;
			best_score = d[best];
		}
		const int z = ::get_diag(i, j, d2, std::max(begin*block_len, j0 - j), std::max(last*block_len, j0 - j), std::min((best + 1)*block_len, j1 - j), j0-j, diags, block_len, log, cutoff, best_score, query_bc);
		if (z < std::numeric_limits<int>::max()) {
			assert(diags.back().len > 0);
			assert(diags.back().j >= 0 && diags.back().subject_end() <= slen);
			assert(diags.back().i >= 0 && diags.back().query_end() <= qlen);
		}
	}
}

void Diag_scores::scan_diags(int d_begin, int d_end, sequence query, sequence subject, const Long_score_profile &qp, const Bias_correction &query_bc, bool log, vector<Diagonal_node> &diags, bool fast)
{
	assert(d_end > d_begin);
	qlen = (int)query.length();
	slen = (int)subject.length();
	const int band = d_end - d_begin;
	this->fast = fast;
	this->d_begin = d_begin;
	this->d_end = d_end;
	i_begin = std::max(0, d_end - 1) - band + 1;
	j_begin = i_begin - d_begin;
	const int j1 = std::min(qlen - d_begin, slen);
	sv_max.clear();
	sv_max.resize(16);
	assert(j1 > j_begin);
	score_buf.init(16, j1 - j_begin);
	local_max.init(16, (j1 - j_begin + block_len - 1) / block_len);

	for (int i = i_begin; i < i_begin + band; i += 16) {

		memset(sv_max.data(), 0, sv_max.size());
		
#ifdef __SSE2__
		scan_cols(qp, subject, i, j_begin, j1, sv_max, log, score_buf, local_max, block_len);
#else
		scan_cols(query, subject, i, j_begin, j1, sv_max, log, score_buf, local_max, block_len);
#endif

		for (int o = 0; o < 16; ++o)
			if (sv_max[o] >= Diag_scores::min_diag_score) {
				if (sv_max[o] >= 255 - score_matrix.bias()) {
					const int s = std::min(i + o, 0), i0 = i + o - s, j0 = j_begin - s;
					score_diagonal2(&query[i0], query_bc, &subject[j0], std::min((int)query.length() - i0, (int)subject.length() - j0), i0, j0, diags, fast ? min_diag_score : min_low_score);
				}
				else
					get_diag(i + o, j_begin, o, j_begin, j1, diags, min_diag_score, log, query_bc);
			}

	}
}