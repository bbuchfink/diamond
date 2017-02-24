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

#include "dp.h"

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
			//cout << 'l' << max[0] << max[1] << max[2] << max[3] << endl;
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

void get_zero_index(const Band::Iterator &d, int begin, int end, int max_score, int &z0, int& z1)
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

void get_diag(int i, int j, const Band::Iterator &d, int begin, int max_score, int zero, vector<Diagonal_node> &diags, int block_len,
	int &second_best_begin, int second_best, const sequence &query, const sequence &subject)
{
	int p = begin, p_end = p + block_len;
	for (; p < p_end; ++p)
		if (d[p] == max_score)
			break;
	diags.push_back(Diagonal_segment(i + zero + 1, j + zero + 1, p - zero, max_score));

	if (second_best_begin >= 0) {
		second_best_begin -= 1;
		const Diagonal_segment second = score_diagonal(&query[i + second_best_begin * 16],
			&subject[j + second_best_begin * 16],
			(second_best - second_best_begin + 1) * 16,
			i + second_best_begin * 16,
			j + second_best_begin * 16);
		if (second.i > diags.back().subject_end() && second.score >= config.min_diag_raw_score)
			diags.push_back(second);
		second_best_begin = -1;
	}
}

void Diag_scores::get_diag(int i, int j, int o, const sequence &query, const sequence &subject, vector<Diagonal_node> &diags)
{
	Band::Iterator d(local_max.diag(o)), d2(score_buf.diag(o));
	int p = 0, p_end = local_max.cols(), best = -1, best_score = -1, second_best_score = -1, second_best = -1, second_best_begin = -1, last_zero = -1, last_block = -1, max_score;
	for (; p < p_end; ++p)
		if ((max_score = d[p]) >= config.min_diag_raw_score) {
			if (last_block == -1) {
				int z0, z1;
				get_zero_index(d2, 0, (p + 1)*block_len, max_score, z0, z1);
				if (z1 > 0) {
					::get_diag(i, j, d2, p * block_len, max_score, z0, diags, block_len, second_best_begin, -1, query, subject);
					max_score = -1;
					last_zero = z1;
				}
				else
					last_zero = z0;
			}
			else {
				int z0, z1;
				get_zero_index(d2, (last_block + 1)*block_len, (p + 1)*block_len, max_score, z0, z1);
				if (z0 > last_zero) {
					::get_diag(i, j, d2, best * block_len, best_score, last_zero, diags, block_len, second_best_begin, second_best, query, subject);
					best = p;
					best_score = max_score;
					last_zero = z0;
				}
				if (z1 > 0) {
					if (max_score > best_score) {
						best = p;
						best_score = max_score;
					}
					::get_diag(i, j, d2, best * block_len, best_score, last_zero, diags, block_len, second_best_begin, second_best, query, subject);
					max_score = -1;
					best = -1;
					best_score = -1;
					last_zero = z1;
				}
			}
			last_block = p;

			if (max_score > best_score) {
				best = p;
				best_score = max_score;
			}
			else if (max_score > d[p - 1]) {
				if (second_best_begin == -1) {
					second_best_begin = p;
					second_best = p;
					second_best_score = max_score;
				}
				else {
					if (max_score > second_best_score) {
						second_best_score = max_score;
						second_best = p;
					}
				}
			}
		}
	if (best >= 0)
		::get_diag(i, j, d2, best * 16, best_score, last_zero, diags, block_len, second_best_begin, second_best, query, subject);
}

void Diag_scores::scan_diags(const Diagonal_segment &diag, sequence query, sequence subject, const Long_score_profile &qp, bool log, vector<Diagonal_node> &diags)
{
	const int d = diag.diag() - band / 2,
		d1 = d + band - 1,
		i = std::max(0, d1) - band + 1,
		j = i - d,
		j1 = std::min((int)query.length() - d, (int)subject.length());
	uint8_t sv_max[band];
	memset(sv_max, 0, band);
	const size_t cells = band * (j1 - j);
	score_buf.init(band, j1 - j);
	//memset(score_buf.data(), 0, score_buf.size());
	local_max.init(band, (j1 - j + block_len - 1) / block_len);
	//memset(local_max.data(), 0, local_max.size());
	scan_cols(qp, subject, i, j, j1, sv_max, log, score_buf.data(), local_max.data(), block_len);
	for (int o = 0; o < band; ++o)
		if (sv_max[o] >= config.min_diag_raw_score) {
			//get_diag(i + o, j, sv_max[o], o);
			if (sv_max[o] >= 255 - score_matrix.bias()) {
				const int s = std::min(i + o, 0);
				diags.push_back(score_diagonal(&query[i + o - s], &subject[j - s], i + o - s, j - s));
			}
			else
				get_diag(i + o, j, o, query, subject, diags);
		}
}