/****
Copyright © 2013-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-2-Clause

#pragma once
#include <array>
#include <ostream>
#include <vector>
#include <string>
#include <assert.h>
#include "basic/value.h"
#include "util/text_buffer.h"
#include "translated_position.h"
#include "util/geo/interval.h"

struct Sequence
{
	static constexpr Letter DELIMITER = DELIMITER_LETTER;
	struct Reversed {};
	struct Hardmasked {};
	Sequence():
		len_ (0),
		data_ (nullptr)
	{ }
	Sequence(const Letter *data, Loc len):
		len_ (len),
		data_ (data)
	{ }
	Sequence(const Letter* data, int64_t len) :
		len_((Loc)len),
		data_(data)
	{ }
	Sequence(const Letter* data, size_t len) :
		len_((Loc)len),
		data_(data)
	{ }
	Sequence(const Letter *begin, const Letter *end) :
		len_(Loc(end - begin)),
		data_(begin)
	{}
	Sequence(const std::vector<Letter> &data):
		len_((Loc)data.size()),
		data_(data.data())
	{}
	Sequence(const Sequence &seq, int from, int to):
		len_(to-from+1),
		data_(seq.data() + from)
	{}
	Loc length() const
	{
		return len_;
	}
	const Letter* data() const
	{
		return data_;
	}
	const Letter* end() const
	{
		return data_ + len_;
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
	Sequence operator+(int d) const {
		return Sequence(data_ + d, len_ - d);
	}
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
	std::string to_string() const {
		std::string s;
		s.resize(len_);
		print(&s[0], 0, len_);
		return s;
	}
	TextBuffer& print(TextBuffer &buf, size_t begin, size_t end, const ValueTraits& value_traits) const
	{
		for (size_t i = begin; i < end; ++i)
#ifdef SEQ_MASK
			buf << value_traits.alphabet[(long)(data_[i] & LETTER_MASK)];
#else
			buf << value_traits.alphabet[(long)data_[i]];
#endif
		return buf;
	}
	std::ostream& print(std::ostream &os, const ValueTraits&v) const
	{
		for (Loc i = 0; i < len_; ++i) {
			long l = (long)data_[i];
			if ((l & 128) == 0)
				os << v.alphabet[l];
			else
				os << (char)tolower(v.alphabet[l & 127]);
		}
		return os;
	}
	TextBuffer& print(TextBuffer &os, const ValueTraits&v) const
	{
		for (Loc i = 0; i < len_; ++i) {
			long l = (long)data_[i];
			if ((l & 128) == 0)
				os << v.alphabet[l];
			else
				os << (char)tolower(v.alphabet[l & 127]);
		}
		return os;
	}
	TextBuffer& print(TextBuffer &os, const ValueTraits&v, Hardmasked) const
	{
		for (Loc i = 0; i < len_; ++i) {
			long l = (long)data_[i];
			if ((l & 128) == 0)
				os << v.alphabet[l];
			else
				os << v.alphabet[(long)v.mask_char];
		}
		return os;
	}
	TextBuffer& print(TextBuffer &os, const ValueTraits&v, Reversed) const
	{
		for (int i = (int)len_ - 1; i >= 0; --i)
			os << v.alphabet[(long)(data_[i] & 127)];
		return os;
	}
	Sequence subseq(int begin, int end) const
	{
		return Sequence(*this, begin, end - 1);
	}
	Sequence subseq(int begin) const {
		return Sequence(data_ + begin, len_ - begin);
	}
	Sequence subseq_clipped(int begin, int end) const {
		return Sequence(*this, std::max(begin, 0), std::min(end, len_));
	}
	friend TextBuffer& operator<<(TextBuffer &buf, const Sequence &s)
	{
		return s.print(buf, value_traits);
	}
	friend std::ostream& operator<<(std::ostream &buf, const Sequence &s)
	{
		return s.print(buf, value_traits);
	}
	static Sequence get_window(const Letter *s, int window)
	{
		const Letter *p = s;
		int n = 0;
		while (*p != Sequence::DELIMITER && n < window) {
			--p;
			++n;
		}
		n = 0;
		while (*s != Sequence::DELIMITER && n < window) {
			++s;
			++n;
		}
		return Sequence(p + 1, Loc(s - p - 1));
	}
	std::vector<Letter> copy() const {
		return std::vector<Letter>(data_, data_ + len_);
	}
	std::vector<Letter> reverse() const;
	void mask(const Interval &i) {
		for (int j = i.begin_; j < i.end_; ++j)
			((Letter*)data_)[j] = value_traits.mask_char;
	}
	bool operator==(const Sequence& s) const {
		if (len_ != s.len_)
			return false;
		for (Loc i = 0; i < len_; ++i)
			if (letter_mask(data_[i]) != letter_mask(s.data_[i]))
				return false;
		return true;
	}
	Loc masked_letters() const {
		Loc n = 0;
		for (Loc i = 0; i < len_; ++i)
			if (letter_mask(data_[i]) == MASK_LETTER)
				++n;
		return n;
	}
	double masked_letter_ratio() const {
		return (double)masked_letters() / len_;
	}
	double length_ratio(const Sequence& seq) const {
		return len_ < seq.len_ ? (double)len_ / seq.len_ : (double)seq.len_ / len_;
	}
	static std::vector<Letter> from_string(const char* str, const ValueTraits&vt = value_traits, int64_t line = 0);

	Loc len_;
	const Letter *data_;
};

struct TranslatedSequence
{

	TranslatedSequence()
	{}

	explicit TranslatedSequence(const Sequence &s1):
		source_(s1)
	{
		translated_[0] = s1;
	}

	TranslatedSequence(const Sequence&source, const Sequence&s1, const Sequence&s2, const Sequence&s3, const Sequence&s4, const Sequence&s5, const Sequence&s6):
		source_(source)
	{
		translated_[0] = s1;
		translated_[1] = s2;
		translated_[2] = s3;
		translated_[3] = s4;
		translated_[4] = s5;
		translated_[5] = s6;
	}

	TranslatedSequence(const Sequence& source, const std::array<std::vector<Letter>, 6>& v):
		source_(source)
	{
		for (int i = 0; i < 6; ++i)
			translated_[i] = Sequence(v.at(i));
	}

	const Sequence& operator[](Frame frame) const
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

	const Sequence& index(unsigned frame) const
	{
		return translated_[frame];
	}

	const Sequence&source() const
	{
		return source_;
	}

	bool in_bounds(const TranslatedPosition &i) const
	{
		return i >= 0 && i < int((*this)[i.frame].length());
	}

	void get_strand(Strand strand, Sequence*dst) const
	{
		int i = strand == FORWARD ? 0 : 3;
		dst[0] = translated_[i++];
		dst[1] = translated_[i++];
		dst[2] = translated_[i];
	}

private:
	
	Sequence source_;
	Sequence translated_[6];

};