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

#ifndef SCORE_PROFILE_H_
#define SCORE_PROFILE_H_

#include <vector>
#include "../basic/sequence.h"
#include "score_vector.h"

using std::vector;

#ifdef __SSE2__

struct sequence_stream
{
	sequence_stream():
		next (buffer_len),
		mask (0)
	{ }
	inline void reset()
	{
		next = buffer_len;
		mask = 0;
	}
	template<typename _score>
	inline const __m128i& get(const typename vector<sequence>::const_iterator &begin,
					   const typename vector<sequence>::const_iterator &end,
					   unsigned pos,
					   const _score&)
	{
		if(next == buffer_len)
			fill<_score>(begin, end, pos);
		return data_[next++];
	}
	template<typename _score>
	inline void fill(const typename vector<sequence>::const_iterator &begin,
		  	  const typename vector<sequence>::const_iterator &end,
		 	  unsigned pos)
	{
		memset(data_, value_traits.mask_char, buffer_len*16);
		unsigned n = 0;
		typename vector<sequence>::const_iterator it (begin);
		assert(pos < it->length());
		const unsigned read_len (std::min(unsigned(buffer_len), static_cast<unsigned>(it->length())-pos));
		while(it < end) {
			const uint8_t *src (reinterpret_cast<const uint8_t*>(it->data()) + pos);
			_score *dest (reinterpret_cast<_score*>(data_) + n);
			int clip (int(pos) - it->clipping_offset_);
			if((mask & (1 << n)) == 0) {
				if(copy_char(src, dest, mask, n, clip))
				if(read_len > 1 && copy_char(src, dest, mask, n, clip))
				if(read_len > 2 && copy_char(src, dest, mask, n, clip))
				if(read_len > 3) copy_char(src, dest, mask, n, clip);
			}
			++it;
			++n;
		}
		next = 0;
	}
	template<typename _score>
	static inline bool copy_char(const uint8_t*& src, _score*& dest, unsigned &mask, unsigned n, int &clip)
	{
		if(clip++ < 0) {
			dest += 16/sizeof(_score);
			++src;
			return true;
		}
		if(*src == 0xff) {
			mask |= 1 << n;
			return false;
		}
		*dest = *(src++) & 0x7f;
		dest += 16/sizeof(_score);
		return true;
	}
	static const unsigned buffer_len = 4;
	__m128i data_[buffer_len];
	unsigned next;
	unsigned mask;
};

template<typename _score>
struct score_profile
{

	inline void set(const __m128i &seq)
	{
		assert(sizeof(data_)/sizeof(score_vector<_score>) >= value_traits.alphabet_size);
		/*unsigned j = 0;
		do {
			data_[j] = score_vector<_score> (j, seq);
			++j;
			data_[j] = score_vector<_score> (j, seq);
			++j;
			data_[j] = score_vector<_score> (j, seq);
			++j;
			data_[j] = score_vector<_score> (j, seq);
			++j;
		} while(j<24);
		data_[j] = score_vector<_score> (j, seq);
		assert(j+1 == Value_traits<_val>::ALPHABET_SIZE);*/
		for (unsigned j = 0; j < value_traits.alphabet_size; ++j)
			data_[j] = score_vector<_score> (j, seq);
	}

	inline const score_vector<_score>& get(Letter i) const
	{
		return data_[(int)i];
	}

	score_vector<_score> data_[25];

};

#endif

struct Long_score_profile
{
	Long_score_profile()
	{}
	Long_score_profile(sequence seq)
	{
		for (unsigned l = 0; l < 25; ++l) {
			const uint8_t *scores = &score_matrix.matrix8u()[l << 5];
			data[l].reserve(seq.length() + 2*padding);
			data[l].insert(data[l].end(), padding, 0);
			for (unsigned i = 0; i < seq.length(); ++i)
				data[l].push_back(scores[(int)seq[i]]);
			data[l].insert(data[l].end(), padding, 0);
		}
	}
	size_t length() const
	{
		return data[0].size() - 2 * padding;
	}
	const uint8_t* get(Letter l, int i) const
	{
		return &data[(int)l][i + padding];
	}
	vector<uint8_t> data[25];
	enum { padding = 32 };
};

#endif /* SCORE_PROFILE_H_ */
