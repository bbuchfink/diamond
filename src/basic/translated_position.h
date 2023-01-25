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

#pragma once
#include <algorithm>
#include "../util/geo/interval.h"
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
		if (!align_mode.query_translated)
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
