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

int Diag_scores::min_diag_score = 22;
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

void scan_cols(const Long_score_profile &qp, sequence s, int i, int i_end, int j, int j_end, uint8_t *sv_max, bool log, uint8_t *buf, uint8_t *local_max, int block_len)
{
#ifdef __SSE2__
	typedef score_vector<uint8_t> Sv;
	const Sv vbias(score_matrix.bias());
	for (int i2 = i; i2 < i_end; i2 += 16) {
		int j2 = std::max(-(i2 - j + 15), j), i3 = j2 + i2 - j;
		uint8_t *local_max_ptr = local_max + (i2 - i) + (j2 - j) / 16 * (i_end - i),
			*buf_ptr = buf + (i2 - i) + (j2 - j)*(i_end - i);
		Sv v, max, global_max;
		for (; j2 < j_end; ++j2, ++i3) {
			assert(j2 >= 0);
			const uint8_t *q = qp.get(s[j2], i3);
			v = v + score_vector<uint8_t>(q);
			v -= vbias;
			max.max(v);
			v.store(buf_ptr);
			buf_ptr += i_end - i;
			if (((j2 - j) & 15) == 15) {
				global_max.max(max);
				max.store(local_max_ptr);
				local_max_ptr += (i_end - i);
				max = score_vector<uint8_t>();
			}
		}
		if (((j2 - j) & 15) != 0) {
			global_max.max(max);
			max.store(local_max_ptr);
		}
		global_max.store(&sv_max[i2 - i]);
	}
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

void Diag_scores::get_diag(int i, int j, int o, int j_begin, int j_end, vector<Diagonal_node> &diags, int cutoff, bool log)
{
	Band::Iterator d(local_max.diag(o)), d2(score_buf.diag(o));
	int p = (j_begin - j) / block_len, p_begin = p, p_end = (j_end - j + block_len - 1) / block_len, best = -1, best_score = -1, last_zero = -1, last_block = -1, max_score, p2_end = score_buf.cols();
	for (; p < p_end; ++p)
		if ((max_score = d[p]) >= cutoff) {
			if (last_block == -1) {
				int z0, z1;
				best = p;
				best_score = max_score;
				if (j_begin == j)
					get_zero_index(d2, 0, std::min((p + 1)*block_len, p2_end), max_score, z0, z1);
				else
					get_zero_index(d2, 0, std::min((p + 1)*block_len, p2_end), max_score, z0, z1, best_score, best);
				if (z1 > 0) {
					const int z = ::get_diag(i, j, d2, best * block_len, best_score, z0, diags, block_len, score_buf.cols(), false, log);
					set_zero(d, d2, z0, z);
					max_score = -1;
					best = -1;
					best_score = -1;
					last_zero = z1;
				}
				else
					last_zero = z0;
			}
			else {
				int z0, z1;
				get_zero_index(d2, (last_block + 1)*block_len, std::min((p + 1)*block_len, p2_end), max_score, z0, z1);
				if (z0 > last_zero) {
					if (best >= 0) {
						const int z = ::get_diag(i, j, d2, best * block_len, best_score, last_zero, diags, block_len, score_buf.cols(), false, log);
						set_zero(d, d2, last_zero, z);
					}
					best = p;
					best_score = max_score;
					last_zero = z0;
				}
				if (z1 > 0) {
					if (max_score > best_score) {
						best = p;
						best_score = max_score;
					}
					const int z = ::get_diag(i, j, d2, best * block_len, best_score, last_zero, diags, block_len, score_buf.cols(), false, log);
					set_zero(d, d2, last_zero, z);
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
		}
	if (best >= 0) {
		const int z = ::get_diag(i, j, d2, best * block_len, best_score, last_zero, diags, block_len, score_buf.cols(), true, log), blockz = z / block_len;
		set_zero(d, d2, last_zero, z);
	}
}

void Diag_scores::scan_ends(unsigned d_idx, vector<Diagonal_node> &diags, bool log, uint8_t *sv_max)
{
	const int dd = diags[d_idx].diag();
	int min_score;
	for (int shift = 1; (min_score = config.gap_open + shift*config.gap_extend + 1) < min_diag_score; shift++)
	{
		if (dd + shift < d_end) {
			const int o = dd + shift - d_begin;
			/*if(sv_max[o] >= min_score)
				get_diag(i_begin+o,j_begin, o, )*/
		}
		if (dd - shift >= d_begin) {

		}
	}
}

void Diag_scores::scan_vicinity(unsigned d_idx, unsigned e_idx, vector<Diagonal_node> &diags, bool log, uint8_t *sv_max)
{
	static const int reverse_diags = 2;
	const int dd = diags[d_idx].diag(), de = diags[e_idx].diag(), jd = diags[d_idx].j, je = diags[e_idx].j, je1 = diags[e_idx].subject_end(), jd1 = diags[d_idx].subject_end(), ld = diags[d_idx].len, le = diags[e_idx].len;
	const int shift = dd - de;
	if (shift >= 0) {
		const int j0 = std::min(je1, jd),
			j1 = std::min(std::max(jd, je1), jd1);
		for (int diag = de + 1; diag < dd; ++diag) {
			const int o = diag - d_begin;
			assert(j0 >= 0);
			assert(j1 <= jd1);
			if(sv_max[o] >= min_low_score)
				get_diag(i_begin + o, j_begin, o, j0, j1, diags, min_low_score, log);
		}
		for (int diag = std::max(de - reverse_diags,d_begin); diag < de; ++diag) {
			const int o = diag - d_begin;
			assert(j0 >= 0);
			assert(j1 <= jd1);
			if (sv_max[o] >= min_low_score)
				get_diag(i_begin + o, j_begin, o, j0, j1, diags, min_low_score, log);
		}
	} else {
		const int jde = jd + shift,
			jde1 = jd1 + shift;
		int j0, l;
		if (jde > je1) {
			j0 = je1;
			l = jde - je1;
		}
		else if (jde > je) {
			j0 = jde;
			l = std::min(je1 - jde, ld);
		}
		else {
			j0 = je;
			l = std::min(jde1 - je, le);
			if (l <= 0)
				return;
		}

		int j = j0;
		for (int diag = de - 1; diag > dd; --diag) {
			const int o = diag - d_begin;
			assert(l >= 0);
			assert(j >= 0);
			assert(j + l <= std::max(diags[d_idx].subject_end(), diags[e_idx].subject_end()));			
			if (sv_max[o] >= min_low_score)
				get_diag(i_begin + o, j_begin, o, j, j + l, diags, min_low_score, log);
			++j;
		}
		for (int diag = de + 1; diag < std::min(de + reverse_diags + 1, d_end); ++diag) {
			const int o = diag - d_begin;
			assert(l >= 0);
			assert(j0 >= 0);
			assert(j0 + l <= std::max(diags[d_idx].subject_end(), diags[e_idx].subject_end()));
			if (sv_max[o] >= min_low_score)
				get_diag(i_begin + o, j_begin, o, j0, j0 + l, diags, min_low_score, log);
		}
	}
}

void Diag_scores::scan_diags(const Diagonal_segment &diag, sequence query, sequence subject, const Long_score_profile &qp, bool log, vector<Diagonal_node> &diags, multimap<int, unsigned> &window)
{
	static const int band = 128, max_dist = 60;
	d_begin = diag.diag() - band / 2;
	d_end = d_begin + band;
	i_begin = std::max(0, d_end - 1) - band + 1;
	j_begin = i_begin - d_begin;
	const int j1 = std::min((int)query.length() - d_begin, (int)subject.length());
	uint8_t sv_max[band];
	memset(sv_max, 0, band);
	score_buf.init(band, j1 - j_begin, true);
	local_max.init(band, (j1 - j_begin + block_len - 1) / block_len, true);
	//scan_cols(qp, subject, i, j, j1, sv_max, log, score_buf.data(), local_max.data(), block_len);
	scan_cols(qp, subject, i_begin, i_begin + band, j_begin, j1, sv_max, log, score_buf.data(), local_max.data(), block_len);
	for (int o = 0; o < band; ++o)
		if (sv_max[o] >= Diag_scores::min_diag_score) {
			if (sv_max[o] >= 255 - score_matrix.bias()) {
				const int s = std::min(i_begin + o, 0), i0 = i_begin + o - s, j0 = j_begin - s;
				score_diagonal(&query[i0], &subject[j0], std::min(query.length() - i0, subject.length() - j0), i0, j0, diags);
				::set_zero(local_max.diag(o), 0, local_max.cols());
			}
			else
				get_diag(i_begin + o, j_begin, o, j_begin, j1, diags, min_diag_score, log);
		}

	std::sort(diags.begin(), diags.end(), Diagonal_segment::cmp_subject);

	const unsigned nodes = diags.size();
	for (unsigned node = 0; node < nodes; ++node) {
		if (diags[node].score < min_diag_score)
			continue;
		if (log)
			cout << "Node n=" << node << endl;
		while (!window.empty() && diags[node].j - window.begin()->first > max_dist)
			window.erase(window.begin());
		for (multimap<int, unsigned>::iterator i = window.begin(); i != window.end(); ++i) {
			if (abs(diags[node].diag() - diags[i->second].diag()) > max_dist)
				continue;
			if (log)
				cout << "Link n=" << i->second << endl;
			scan_vicinity(node, i->second, diags, log, sv_max);
			/*if (e.subject_end() - (d.subject_end() - std::min(e.diag() - d.diag(), 0)) >= 11) {

			}*/
		}
		window.insert(std::make_pair(diags[node].subject_end(), node));
	}

	window.clear();
}