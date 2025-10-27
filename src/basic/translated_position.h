/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "util/geo/interval.h"
#include "value.h"

enum Strand { FORWARD = 0, REVERSE = 1 };

struct Frame
{
	Frame()
	{}
	Frame(Strand strand, int offset):
		offset(offset),
		strand(strand)		
	{}
	explicit Frame(int index) :
		offset(index % 3),
		strand(index < 3 ? FORWARD : REVERSE)		
	{}
	int index() const
	{
		return strand * 3 + offset;
	}
	int signed_frame() const
	{
		return (offset + 1) * (strand == FORWARD ? 1 : -1);
	}
	int length(int dna_len) const
	{
		return std::max((dna_len - offset) / 3, 0);
	}
	int offset;
	Strand strand;
};

struct TranslatedPosition
{

	TranslatedPosition()
	{}

	TranslatedPosition(int translated, Frame frame):
		frame(frame),
		translated(translated)		
	{
	}
	
	TranslatedPosition(int in_strand, Strand strand):
		frame(strand, in_strand % 3),
		translated(in_strand_to_translated(in_strand))
	{}

	operator int() const
	{
		return translated;
	}

	TranslatedPosition& operator++()
	{
		++translated;
		return *this;
	}

	TranslatedPosition operator+(int x) const
	{
		return TranslatedPosition(translated + x, frame);
	}

	TranslatedPosition operator-(int x) const
	{
		return TranslatedPosition(translated - x, frame);
	}

	void shift_forward()
	{
		++frame.offset;
		if (frame.offset == 3) {
			frame.offset = 0;
			++translated;
		}
	}

	void shift_back()
	{
		--frame.offset;
		if (frame.offset == -1) {
			frame.offset = 2;
			--translated;
		}
	}

	void shift_forward(int k) {
		while (k-- > 0)
			shift_forward();
	}

	int frame_shift(const TranslatedPosition &x) const
	{
		static const int frameshift[3][3] = { { 0, 1, -1 },{ -1, 0, 1 },{ 1, -1, 0 } };
		return frameshift[frame.offset][x.frame.offset];
	}

	int absolute(int dna_len) const
	{
        if(!frame.offset &&  align_mode.mode == AlignMode::blastn)
           return dna_len - 1 - translated;
		if (!align_mode.query_translated && frame.strand == FORWARD)
			return translated;
		return oriented_position(in_strand(), frame.strand, dna_len);
	}

	static Interval absolute_interval(const TranslatedPosition &begin, const TranslatedPosition &end, int dna_len)
	{
		if (begin.frame.strand == FORWARD)
			return Interval(begin.in_strand(), end.in_strand());
		else
			return Interval(oriented_position(end.in_strand() - 1, REVERSE, dna_len), oriented_position(begin.in_strand() - 1, REVERSE, dna_len));
	}

	static int in_strand_to_translated(int in_strand)
	{
		if (align_mode.query_translated)
			return in_strand / 3;
		else
			return in_strand;
	}

	static int translated_to_in_strand(int translated, Frame frame)
	{
		if (align_mode.query_translated)
			return frame.offset + 3 * translated;
		else
			return translated;
	}

	int in_strand() const
	{
		return translated_to_in_strand(translated, frame);
	}

	static int oriented_position(int pos, Strand strand, int dna_len)
	{
		return strand == FORWARD ? pos : dna_len - pos - 1;
	}

	static int absolute_to_translated(int src, Frame frame, int dna_len, bool translated)
	{
		if (!translated)
			return src;
		return in_strand_to_translated(oriented_position(src, frame.strand, dna_len));
	}

	static int translated_to_absolute(int translated, Frame frame, int dna_len)
	{
		return oriented_position(translated_to_in_strand(translated, frame), frame.strand, dna_len);
	}

	friend std::ostream& operator<<(std::ostream &s, const TranslatedPosition &a)
	{
		s << a.frame.offset << ' ' << a.translated;
		return s;
	}

	Frame frame;
	int translated;

};
