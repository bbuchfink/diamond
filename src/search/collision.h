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

#ifndef COLLISION_H_
#define COLLISION_H_

template<typename _val>
inline unsigned letter_match(const _val query, const _val subject)
{
	if(query != Value_traits<_val>::MASK_CHAR && Reduction<_val>::reduction(query) == Reduction<_val>::reduction(mask_critical(subject)))
		return 1;
	else
		return 0;
}

inline bool match_shape_mask(const uint64_t mask, const uint64_t shape_mask)
{
	return (mask & shape_mask) == shape_mask;
}

template<typename _val>
inline bool is_lower_chunk(const _val *subject, unsigned sid)
{
	uint64_t seed;
	shape_config::get().get_shape(sid).set_seed(seed, subject);
	return current_range.lower(seed_partition(seed));
}

template<typename _val>
inline bool is_lower_or_equal_chunk(const _val *subject, unsigned sid)
{
	uint64_t seed;
	shape_config::get().get_shape(sid).set_seed(seed, subject);
	return current_range.lower_or_equal(seed_partition(seed));
}

inline bool need_lookup(unsigned sid)
{
	//return sid != 0 || current_query_chunk != 0;
	return sid != 0;
}

template<typename _val>
inline bool is_low_freq(const _val *subject, unsigned sid)
{ return shape_config::get().get_shape(sid).is_low_freq(subject); }

template <typename _val, typename _pos>
inline bool shape_collision_right(uint64_t mask, uint64_t shape_mask, const _val *subject, unsigned sid)
{
	if(!match_shape_mask(mask, shape_mask)) return false;
	return is_lower_chunk(subject, sid)
			&& is_low_freq(subject, sid)
			&& (!get_critical(*subject) || (need_lookup(sid) && !ref_masking.get(subject, sid)));
}

template <typename _val, typename _pos>
inline bool shape_collision_left(uint64_t mask, uint64_t shape_mask, const _val *subject, unsigned sid, bool chunked)
{
	if(!match_shape_mask(mask, shape_mask)) return false;
	return (!chunked || is_lower_or_equal_chunk(subject, sid))
			&& is_low_freq(subject, sid)
			&& (!get_critical(*subject) || (need_lookup(sid) && !ref_masking.get(subject, sid)));
}

template <typename _val, typename _pos>
inline bool previous_shape_collision(uint64_t mask, uint64_t shape_mask, const _val *subject, unsigned sid)
{
	if(!match_shape_mask(mask, shape_mask)) return false;
	return is_low_freq(subject, sid)
			&& (!get_critical(*subject) || !ref_masking.get(subject, sid));
}

template<typename _val, typename _pos>
bool is_primary_hit2(const _val *query,
					const _val *subject,
					const unsigned seed_offset,
					const unsigned sid,
					const unsigned len,
					Statistics &stat)
{
	assert(len > 0 && len <= program_options::window*2);
	unsigned mask (0);
	const bool chunked (program_options::lowmem > 1);
	unsigned i = len;

	do {
		--i;
		mask <<= 1;
		mask |= letter_match(query[i], subject[i]);
		for(unsigned j=0;j<sid;++j)
			if(previous_shape_collision<_val,_pos>(mask, shape_config::instance.get_shape(j).mask_, &subject[i], j, stat))
				return false;
		if(chunked && i > seed_offset && shape_collision_right<_val,_pos>(mask, shape_config::instance.get_shape(sid).mask_, &subject[i], sid, stat))
			return false;
	} while (i > seed_offset);

	if(i == 0)
		return true;

	do {
		--i;
		mask <<= 1;
		mask |= letter_match(query[i], subject[i]);
		for(unsigned j=0;j<sid;++j)
			if(previous_shape_collision<_val,_pos>(mask, shape_config::instance.get_shape(j).mask_, &subject[i], j, stat))
				return false;
		if(shape_collision_left<_val,_pos>(mask, shape_config::instance.get_shape(sid).mask_, &subject[i], sid, chunked, stat))
			return false;
	} while (i > 0);

	return true;
}

template<typename _val, typename _pos>
bool is_primary_hit(const _val *query,
					const _val *subject,
					const unsigned seed_offset,
					const unsigned sid,
					const unsigned len)
{
	assert(len > 0 && len <= program_options::window*2);
	const bool chunked (program_options::lowmem > 1);
	uint64_t mask = reduced_match32(query, subject, len);
	unsigned i = 0;
	uint64_t current_mask = shape_config::instance.get_shape(sid).mask_;
	unsigned shape_len =  len - shape_config::instance.get_shape(0).length_ + 1;
	while(i < shape_len) {
		if(len-i > 32)
			mask |= reduced_match32(query+32,subject+32,len-i-32) << 32;
		for(unsigned j=0;j<32 && i<shape_len;++j) {
			assert(&subject[j] >= ref_seqs<_val>::data_->data(0) && &subject[j] <= ref_seqs<_val>::data_->data(ref_seqs<_val>::data_->raw_len()-1));
			for(unsigned k=0;k<sid;++k)
				if(previous_shape_collision<_val,_pos>(mask, shape_config::instance.get_shape(k).mask_, &subject[j], k))
					return false;
			if(i < seed_offset && shape_collision_left<_val,_pos>(mask, current_mask, &subject[j], sid, chunked))
				return false;
			if(chunked && i > seed_offset && shape_collision_right<_val,_pos>(mask, current_mask, &subject[j], sid))
				return false;
			++i;
			mask >>= 1;
		}
		query += 32;
		subject += 32;
	}
	return true;
}

#endif /* COLLISION_H_ */
