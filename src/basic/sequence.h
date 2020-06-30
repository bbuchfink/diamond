/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include "../basic/value.h"
#include "../util/binary_buffer.h"
#include "../util/text_buffer.h"
#include "translated_position.h"
#include "../util/interval.h"

using std::vector;

struct sequence
{
	static constexpr Letter DELIMITER = DELIMITER_LETTER;
	struct Reversed {};
	struct Hardmasked {};
	sequence():
		len_ (0),
		clipping_offset_ (0),
		data_ (0)
	{ }
	sequence(const Letter *data, size_t len, int clipping_offset = 0):
		len_ (len),
		clipping_offset_ (clipping_offset),
		data_ (data)
	{ }
	sequence(const Letter *begin, const Letter *end) :
		len_(end - begin),
		clipping_offset_(0),
		data_(begin)
	{}
	sequence(const vector<Letter> &data):
		len_(data.size()),
		clipping_offset_(0),
		data_(data.data())
	{}
	sequence(const sequence &seq, int from, int to):
		len_(to-from+1),
		clipping_offset_(0),
		data_(seq.data() + from)
	{}
	size_t length() const
	{
		return len_;
	}
	size_t clipped_length() const
	{
		return len_ - clipping_offset_;
	}
	size_t aligned_clip(unsigned padding) const
	{
		return clipping_offset_ > (int)padding ? clipping_offset_ - padding : 0;
	}
	const Letter* data() const
	{
		return data_;
	}
	const Letter* end() const
	{
		return data_ + len_;
	}
	const Letter* clipped_data() const
	{
		return data_ + clipping_offset_;
	}
	const Letter* aligned_data(unsigned padding) const
	{
		return data_ + padding;
	}
	Letter operator [](size_t i) const
	{
#ifdef SEQ_MASK
		return data_[i] & LETTER_MASK;
#else
		return data_[i];
#endif
	}
	bool empty() const
	{ return len_ == 0; }
	size_t print(char *ptr, unsigned begin, unsigned len) const
	{
		for(unsigned i=begin;i<begin+len;++i)
#ifdef SEQ_MASK
			* (ptr++) = to_char(data_[i] & LETTER_MASK);
#else
			*(ptr++) = to_char(data_[i]);
#endif
		return len;
	}
	TextBuffer& print(TextBuffer &buf, size_t begin, size_t end, const Value_traits& value_traits) const
	{
		for (size_t i = begin; i < end; ++i)
#ifdef SEQ_MASK
			buf << value_traits.alphabet[(long)(data_[i] & LETTER_MASK)];
#else
			buf << value_traits.alphabet[(long)data_[i]];
#endif
		return buf;
	}
	std::ostream& print(std::ostream &os, const Value_traits &v) const
	{
		for (unsigned i = 0; i < len_; ++i) {
			long l = (long)data_[i];
			if ((l & 128) == 0)
				os << v.alphabet[l];
			else
				os << (char)tolower(v.alphabet[l & 127]);
		}
		return os;
	}
	TextBuffer& print(TextBuffer &os, const Value_traits &v) const
	{
		for (unsigned i = 0; i < len_; ++i) {
			long l = (long)data_[i];
			if ((l & 128) == 0)
				os << v.alphabet[l];
			else
				os << (char)tolower(v.alphabet[l & 127]);
		}
		return os;
	}
	TextBuffer& print(TextBuffer &os, const Value_traits &v, Hardmasked) const
	{
		for (unsigned i = 0; i < len_; ++i) {
			long l = (long)data_[i];
			if ((l & 128) == 0)
				os << v.alphabet[l];
			else
				os << v.alphabet[(long)v.mask_char];
		}
		return os;
	}
	TextBuffer& print(TextBuffer &os, const Value_traits &v, Reversed) const
	{
		for (int i = (int)len_ - 1; i >= 0; --i)
			os << v.alphabet[(long)(data_[i] & 127)];
		return os;
	}
	sequence subseq(int begin, int end) const
	{
		return sequence(*this, begin, end - 1);
	}
	friend TextBuffer& operator<<(TextBuffer &buf, const sequence &s)
	{
		return s.print(buf, value_traits);
	}
	friend std::ostream& operator<<(std::ostream &buf, const sequence &s)
	{
		return s.print(buf, value_traits);
	}
	static sequence get_window(const Letter *s, int window)
	{
		const Letter *p = s;
		int n = 0;
		while (*p != sequence::DELIMITER && n < window) {
			--p;
			++n;
		}
		n = 0;
		while (*s != sequence::DELIMITER && n < window) {
			++s;
			++n;
		}
		return sequence(p + 1, s - p - 1);
	}
	std::vector<Letter> copy() const {
		std::vector<Letter> v;
		v.reserve(len_);
		for (size_t i = 0; i < len_; ++i)
			v.push_back(data_[i]);
		return v;
	}
	void mask(const interval &i) {
		for (int j = i.begin_; j < i.end_; ++j)
			((Letter*)data_)[j] = value_traits.mask_char;
	}
	/*friend std::ostream& operator<<(std::ostream &os, const sequence &s)
	{
		std::cout << "co = " << s.clipping_offset_ << std::endl;
		for(unsigned i=s.clipping_offset_;i<s.len_;++i) {
			if(s.data_[i] == 24)
				break;
			os << mask_critical(s.data_[i]);
		}
		return os;
	}*/
	static vector<Letter> from_string(const char* str, const Value_traits &vt = value_traits);
	size_t			len_;
	int				clipping_offset_;
	const Letter *data_;
};

struct TranslatedSequence
{

	TranslatedSequence()
	{}

	explicit TranslatedSequence(const sequence &s1):
		source_(s1)
	{
		translated_[0] = s1;
	}

	TranslatedSequence(const sequence &source, const sequence &s1, const sequence &s2, const sequence &s3, const sequence &s4, const sequence &s5, const sequence &s6):
		source_(source)
	{
		translated_[0] = s1;
		translated_[1] = s2;
		translated_[2] = s3;
		translated_[3] = s4;
		translated_[4] = s5;
		translated_[5] = s6;
	}

	TranslatedSequence(const sequence &source, const vector<Letter> v[6]):
		source_(source)
	{
		for (int i = 0; i < 6; ++i)
			translated_[i] = sequence(v[i]);
	}

	const sequence& operator[](Frame frame) const
	{
		return translated_[frame.index()];
	}

	Letter operator[](const TranslatedPosition &i) const
	{
		return (*this)[i.frame][i];
	}

	Letter operator()(int in_strand, Strand strand) const
	{
		assert(in_strand < (int)source_.length() - 2);
		return translated_[in_strand % 3 + (strand == FORWARD ? 0 : 3)][in_strand / 3];
	}

	const sequence& index(unsigned frame) const
	{
		return translated_[frame];
	}

	const sequence &source() const
	{
		return source_;
	}

	bool in_bounds(const TranslatedPosition &i) const
	{
		return i >= 0 && i < int((*this)[i.frame].length());
	}

	void get_strand(Strand strand, sequence *dst) const
	{
		int i = strand == FORWARD ? 0 : 3;
		dst[0] = translated_[i++];
		dst[1] = translated_[i++];
		dst[2] = translated_[i];
	}

private:
	
	sequence source_;
	sequence translated_[6];

};


#endif /* SEQUENCE_H_ */
