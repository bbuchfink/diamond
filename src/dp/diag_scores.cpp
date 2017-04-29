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

#include <map>
#include "dp.h"

using std::multimap;

int Diag_scores::min_diag_score = 19;
int Diag_scores::min_low_score = 13;

Diagonal_segment score_diagonal(const Letter *query, const Letter *subject, int qbegin, int jbegin)
{
	int i = 0, j = 0, max_score = 0, score = 0, begin = 0, end = 0;
	while (query[i] != '\xff' && subject[i] != '\xff') {
		score += score_matrix(query[i], subject[i]);
		if (score <= 0) {
			score = 0;
			j = i + 1;
		}
		if (score > max_score) {
			max_score = score;
			begin = j;
			end = i + 1;
		}
		++i;
	}
	return Diagonal_segment(qbegin + begin, jbegin + begin, end - begin, max_score);
}

void score_diagonal(const Letter *query, const Letter *subject, int len, int qbegin, int jbegin, vector<Diagonal_node> &diags)
{
	int i = 0, j = 0, max_score = 0, score = 0, begin = 0, end = 0, l = 0;
	while (l<len) {
		score += score_matrix(query[i], subject[i]);
		if (score <= 0) {
			if (max_score >= Diag_scores::min_low_score) {
				diags.push_back(Diagonal_node(qbegin + begin, jbegin + begin, end - begin, max_score));
				max_score = 0;
				score_diagonal(query + end, subject + end, i - end, qbegin + end, jbegin + end, diags);
			}
			score = 0;
			j = i + 1;
		}
		else if (score > max_score) {
			max_score = score;
			begin = j;
			end = i + 1;
		}
		++i;
		++l;
	}
	if (max_score >= Diag_scores::min_low_score) {
		diags.push_back(Diagonal_node(qbegin + begin, jbegin + begin, end - begin, max_score));
		score_diagonal(query + end, subject + end, i - end, qbegin + end, jbegin + end, diags);
	}
}

void score_diagonal2(const Letter *query, const Bias_correction &query_bc, const Letter *subject, int len, int qbegin, int jbegin, vector<Diagonal_node> &diags, int cutoff, size_t &cells)
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

Diagonal_segment score_diagonal(const Letter *query, const Letter *subject, int n, int qbegin, int jbegin)
{
	int i = 0, j = 0, max_score = 0, score = 0, begin = 0, end = 0;
	while (i < n && query[i] != '\xff' && subject[i] != '\xff') {
		score += score_matrix(query[i], subject[i]);
		if (score <= 0) {
			score = 0;
			j = i + 1;
		}
		if (score > max_score) {
			max_score = score;
			begin = j;
			end = i + 1;
		}
		++i;
	}
	return Diagonal_segment(qbegin + begin, jbegin + begin, end - begin, max_score);
}

void set_global_max(score_vector<uint8_t> *max, score_vector<uint8_t> *global_max, uint8_t *&local_max)
{
#ifdef __SSE2__
	global_max[0].max(max[0]);
	max[0].store(local_max);
	max[0] = score_vector<uint8_t>();
	local_max += 16;
	global_max[1].max(max[1]);
	max[1].store(local_max);
	max[1] = score_vector<uint8_t>();
	local_max += 16;
	global_max[2].max(max[2]);
	max[2].store(local_max);
	max[2] = score_vector<uint8_t>();
	local_max += 16;
	global_max[3].max(max[3]);
	max[3].store(local_max);
	max[3] = score_vector<uint8_t>();
	local_max += 16;
#endif
}

void scan_cols(const Long_score_profile &qp, sequence s, int i, int j, int j_end, uint8_t *sv_max, bool log, uint8_t *buf, uint8_t *local_max, int block_len)
{
#ifdef __SSE2__
	typedef score_vector<uint8_t> Sv;
	const Sv vbias(score_matrix.bias());
	Sv v[4], max[4], global_max[4];
	int n = 0;
	for (; j < j_end; ++j) {
		const uint8_t *q = qp.get(s[j], i);
		v[0] = v[0] + score_vector<uint8_t>(q);
		v[0] -= vbias;
		max[0].max(v[0]);
		_mm_storeu_si128((__m128i*)buf, v[0].data_);
		q += 16;
		buf += 16;
		v[1] = v[1] + score_vector<uint8_t>(q);
		v[1] -= vbias;
		max[1].max(v[1]);
		_mm_storeu_si128((__m128i*)buf, v[1].data_);
		q += 16;
		buf += 16;
		v[2] = v[2] + score_vector<uint8_t>(q);
		v[2] -= vbias;
		max[2].max(v[2]);
		_mm_storeu_si128((__m128i*)buf, v[2].data_);
		q += 16;
		buf += 16;
		v[3] = v[3] + score_vector<uint8_t>(q);
		v[3] -= vbias;
		max[3].max(v[3]);
		_mm_storeu_si128((__m128i*)buf, v[3].data_);
		buf += 16;
		//cout << 's' << v[0] << v[1] << v[2] << v[3] << endl;
		if ((n & 15) == 15) {
			//cout << 'l' << max[0] << max[1] << max[2] << max[3] << endl;7
			set_global_max(max, global_max, local_max);
		}
		++i;
		++n;
	}
	if (n % block_len != 0)
		set_global_max(max, global_max, local_max);
	global_max[0].store(sv_max);
	global_max[1].store(sv_max + 16);
	global_max[2].store(sv_max + 32);
	global_max[3].store(sv_max + 48);
	//cout << 'g' << global_max[0] << global_max[1] << global_max[2] << global_max[3] << endl;
#endif
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

#endif
}

void get_zero_index(Band::Iterator &d, int begin, int end, int max_score, int &z0, int& z1)
{
	z0 = z1 = -1;
	for (int i = end - 1; i >= begin; --i) {
		if (d[i] == (uint8_t)0 && z1 == -1)
			z1 = i;
		else if (d[i] == max_score) {
			for (i -= 1; i >= begin; --i) {
				if (d[i] == (uint8_t)0) {
					z0 = i;
					return;
				}
			}
			return;
		}
	}
}

void get_zero_index(Band::Iterator &d, int begin, int end, int max_score, int &z0, int& z1, int &best_score, int &best)
{
	z0 = z1 = -1;
	for (int i = end - 1; i >= begin; --i) {
		if (d[i] == (uint8_t)0 && z1 == -1)
			z1 = i;
		else if (d[i] == max_score) {
			for (i -= 1; i >= begin; --i) {
				if (d[i] > max_score) {
					max_score = d[i];
					best_score = max_score;
					best = i / Diag_scores::block_len;
				}
				else if (d[i] == (uint8_t)0) {
					z0 = i;
					return;
				}
			}
			return;
		}
	}
}

int get_diag(int i, int j, Band::Iterator &d, int begin, int max_score, int zero, vector<Diagonal_node> &diags, int block_len, int cols, bool extend, bool log)
{
	int p = begin, p_end = p + block_len;
	for (; p < p_end; ++p)
		if (d[p] == max_score) {
			if (extend) {
				for (int q = p + 1; q < cols && d[q] > 0;++q)
					if (d[q] > max_score) {
						max_score = d[q];
						p = q;
					}
			}
			break;
		}
	diags.push_back(Diagonal_segment(i + zero + 1, j + zero + 1, p - zero, max_score));
	if (log)
		cout << diags.back() << endl;

	int low = max_score, low_pos = p, high = 0, high_pos = p;
	++p;
	while (p < cols && d[p] > 0) {
		if (d[p] < low) {
			if (high >= Diag_scores::min_diag_score) {
				diags.push_back(Diagonal_segment(i + low_pos + 1, j + low_pos + 1, high_pos - low_pos, high));
				if (log)
					cout << diags.back() << endl;
			}
			high = 0;
			high_pos = p;
			low = d[p];
			low_pos = p;
		}
		if (d[p] - low > high) {
			high = d[p] - low;
			high_pos = p;
		}
		++p;
	}
	if (high >= Diag_scores::min_diag_score) {
		diags.push_back(Diagonal_segment(i + low_pos + 1, j + low_pos + 1, high_pos - low_pos, high));
		if (log)
			cout << diags.back() << endl;
	}
	return p;
}

int get_score_idx(Band::Iterator &d, int begin, int end, int score)
{
	int i;
	for (i = begin; i < end; ++i)
		if (d[i] == score)
			return i;
	return i;
}

int get_max_idx(Band::Iterator &d, int begin, int end)
{
	assert(begin >= 0 && begin < end);
	int i = begin, s = d[begin];
	for (int j=i+1; j < end;++j)
		if (d[j] > s) {
			i = j;
			s = d[j];
		}
	return i;
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

int get_diag(int i, int j, Band::Iterator &d, int begin, int last, int end, int d0, vector<Diagonal_node> &diags, int block_len, bool log, int cutoff, int best_score, size_t &cells, const Bias_correction &query_bc)
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
		if (score >= cutoff || (d[p1]==best_score &&  score > 0)) {
		//if (score >= cutoff) {
			assert(i + p0 + 1 >= 0);
			assert(j + p0 + 1 >= 0);
			Diagonal_segment diag(Diagonal_segment(i + p0 + 1, j + p0 + 1, p1 - p0, score));
			//if (diag.score + query_bc(diag) >= cutoff) {
			diag.score = diag.score; // +query_bc(diag);
				diags.push_back(diag);
				cells += p1 - p0;
				assert(p0 + 1 >= 0);
				z = p0 + 1;
				if (log)
					cout << diags.back() << endl;
			//}
		}
		end = p0;
	}
	assert(z >= last);
	return z;
}

void set_zero(Band::Iterator &i, int begin, int end)
{
	for (int j = begin; j < end; ++j)
		i[j] = (uint8_t)0;
}

bool have_score(Band::Iterator &i, int begin, int end, int score)
{
	for (int k = begin; k < end; ++k)
		if (i[k] == score)
			return true;
	return false;
}

void Diag_scores::set_active(int o, int begin, int end)
{
	if (!fast) {
		vector<bool>::iterator i = active.begin() + o*local_max.cols() + begin, j = i + end - begin;
		for (; i < j; ++i)
			*i = true;
	}
}

bool Diag_scores::is_active(int o, int i) const
{
	return fast ? false : active[o*local_max.cols() + i];
}

void Diag_scores::set_zero(Band::Iterator &d, Band::Iterator d2, int begin, int end)
{
	const int b0 = begin / Diag_scores::block_len, b1 = end / Diag_scores::block_len;
	if (have_score(d2, begin + 1, std::min((b0 + 1)*Diag_scores::block_len, end), d[b0])) {
		assert(b0 >= 0 && b0 < local_max.cols());
		d[b0] = 0;
	}
	if (have_score(d2, b1*Diag_scores::block_len, end, d[b1])) {
		assert(b1 >= 0 && b1 < local_max.cols());
		d[b1] = 0;
	}
	assert(b0 + 1 >= 0 && b1 <= local_max.cols());
	::set_zero(d, b0 + 1, b1);
}

void Diag_scores::get_diag(int i, int j, int o, int j_begin, int j_end, vector<Diagonal_node> &diags, int cutoff, bool log, size_t &cells, const Bias_correction &query_bc)
{
	Band::Iterator d(local_max.diag(o)), d2(score_buf.diag(o));
	const int diag = i - j,
		j0 = dj0(diag),
		j1 = dj1(diag),
		b0 = (j0 - j) / block_len,
		b1 = (j1 - j + block_len - 1) / block_len;
	int p = std::max((j_begin - j) / block_len, b0), p_begin = p, p_end = std::min((j_end - j + block_len - 1) / block_len, b1), best = -1, best_score = -1, begin = -1, max_score, last = p;
	while (last > b0 && !is_active(o, last - 1))
		--last;
	for (; p < p_end; ++p) {
		if (!is_active(o, p)
			&& (max_score = d[p]) >= cutoff
			&& (p == 0 || max_score > d[p - 1])) {
			if (begin == -1)
				begin = p;
			best = p;
			best_score = max_score;
		}
		else if (begin != -1) {
			const int z = ::get_diag(i, j, d2, std::max(begin*block_len, j0 - j), std::max(last*block_len, j0 - j), std::min((best + 1)*block_len, j1 - j), j0-j, diags, block_len, log, cutoff, best_score, cells, query_bc);
			if (z < std::numeric_limits<int>::max()) {
				assert(diags.back().len > 0);
				assert(diags.back().j >= 0 && diags.back().subject_end() <= slen);
				assert(diags.back().i >= 0 && diags.back().query_end() <= qlen);
				set_active(o, z / block_len, best + 1);
				last = best + 1;
			}
			begin = -1;
			best = -1;
		}
		if (is_active(o, p))
			last = p + 1;
	}

	if (begin != -1) {
		if (best == p_end - 1) {
			for (; best<b1 && !is_active(o, best) && d[best] >= cutoff && (best == 0 || d[best] > d[best - 1]); ++best);
			//set_active(o, p_end, best);
			best -= 1;
			best_score = d[best];
		}
		const int z = ::get_diag(i, j, d2, std::max(begin*block_len, j0 - j), std::max(last*block_len, j0 - j), std::min((best + 1)*block_len, j1 - j), j0-j, diags, block_len, log, cutoff, best_score, cells, query_bc);
		if (z < std::numeric_limits<int>::max()) {
			assert(diags.back().len > 0);
			assert(diags.back().j >= 0 && diags.back().subject_end() <= slen);
			assert(diags.back().i >= 0 && diags.back().query_end() <= qlen);
			set_active(o, z / block_len, best + 1);
		}
	}
}

size_t Diag_scores::scan_diags(int d_begin, int d_end, sequence query, sequence subject, const Long_score_profile &qp, const Bias_correction &query_bc, bool log, vector<Diagonal_node> &diags, bool fast)
{
	assert(d_end > d_begin);
	qlen = (int)query.length();
	slen = (int)subject.length();
	d_begin = std::max(d_begin, -((int)subject.length() - 1));
	d_end = std::min(d_end, (int)query.length());
	d_end += (d_end - d_begin) % 16 == 0 ? 0 : 16 - (d_end - d_begin) % 16;
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
	size_t cells = 0;

	for (int i = i_begin; i < i_begin + band; i += 16) {

		memset(sv_max.data(), 0, sv_max.size());
		
		scan_cols(qp, subject, i, j_begin, j1, sv_max, log, score_buf, local_max, block_len);
		//scan_cols(query, subject, i, j_begin, j1, sv_max, log, score_buf, local_max, block_len);

		for (int o = 0; o < 16; ++o)
			if (sv_max[o] >= Diag_scores::min_diag_score) {
				if (sv_max[o] >= 255 - score_matrix.bias()) {
					const int s = std::min(i + o, 0), i0 = i + o - s, j0 = j_begin - s;
					score_diagonal2(&query[i0], query_bc, &subject[j0], std::min((int)query.length() - i0, (int)subject.length() - j0), i0, j0, diags, fast ? min_diag_score : min_low_score, cells);
					//::set_zero(local_max.diag(o), 0, local_max.cols());
					set_active(o, 0, local_max.cols());
				}
				else
					get_diag(i + o, j_begin, o, j_begin, j1, diags, min_diag_score, log, cells, query_bc);
			}

	}
	
	//cout << "dens=" << (double)cells / (double)band / (double)(j1 - j_begin) << endl;

	return cells;

}