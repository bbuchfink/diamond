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
#include "../util/binary_buffer.h"
#include "../basic/value.h"
#include "sequence.h"

typedef enum { op_match = 0, op_insertion = 1, op_deletion = 2, op_substitution = 3, op_frameshift_forward = 4, op_frameshift_reverse = 5 } Edit_operation;

struct Reversed {};

struct Packed_operation
{
	enum { OP_BITS = 2, COUNT_BITS = 8 - OP_BITS, MAX_COUNT = (1 << COUNT_BITS) - 1 };
	Packed_operation(uint8_t code):
		code (code)
	{ }
	Packed_operation(Edit_operation op, unsigned count):
		code ((op<<COUNT_BITS) | count)
	{ }
	Packed_operation(Edit_operation op, Letter v) :
		code((op << COUNT_BITS) | (int)v)
	{ }
	operator uint8_t() const
	{ return code; }
	Edit_operation op() const
	{
		uint8_t o = code >> COUNT_BITS;
		switch (o) {
		case op_substitution:
			switch (letter()) {
			case (Letter)AMINO_ACID_COUNT:
				return op_frameshift_reverse;
			case (Letter)(AMINO_ACID_COUNT + 1) :
				return op_frameshift_forward;
			default:
				return op_substitution;
			}
		default:
			return (Edit_operation)o;
		}
	}
	unsigned count() const
	{
		switch (op()) {
		case op_match:
		case op_insertion:
			return code & MAX_COUNT;
		default:
			return 1;
		}
	}
	Letter letter() const
	{
		return code&MAX_COUNT;
	}
	static Packed_operation terminator()
	{ return Packed_operation(op_match, 0u); }
	static Packed_operation frameshift_forward()
	{
		return Packed_operation(op_substitution, (unsigned)AMINO_ACID_COUNT + 1);
	}
	static Packed_operation frameshift_reverse()
	{
		return Packed_operation(op_substitution, (unsigned)AMINO_ACID_COUNT);
	}
	uint8_t code;
};

struct Combined_operation
{
	Edit_operation op;
	unsigned count;
	Letter letter;
};

struct Packed_transcript
{

	struct Const_iterator
	{
		Const_iterator(const Packed_operation *op):
			ptr_ (op)
		{ gather(); }
		bool good() const
		{ return *ptr_ != Packed_operation::terminator(); }
		Const_iterator& operator++()
		{ ++ptr_; gather(); return *this; }
		const Combined_operation& operator*() const
		{ return op_; }
		const Combined_operation* operator->() const
		{ return &op_; }
	private:
		void gather()
		{
			if(!good())
				return;
			op_.op = ptr_->op();
			if(op_.op == op_deletion || op_.op == op_substitution || op_.op == op_frameshift_forward || op_.op==op_frameshift_reverse) {
				op_.letter = ptr_->letter();
				op_.count = 1;
			} else {
				op_.count = 0;
				do {
					op_.count += ptr_->count();
					++ptr_;
				} while(good() && ptr_->op() == op_.op);
				--ptr_;
			}
		}
		const Packed_operation *ptr_;
		Combined_operation op_;
	};

	void read(BinaryBuffer::Iterator &it)
	{
		data_.clear();
		uint8_t code;
		do {
			it >> code;
			data_.push_back(code);
		} while (code != Packed_operation::terminator());
	}

	Const_iterator begin() const
	{ return Const_iterator (data_.data()); }

	const std::vector<Packed_operation>& data() const
	{ return data_; }

	const Packed_operation* ptr() const
	{
		return &data_[0];
	}

	void push_back(Edit_operation op)
	{
		if (op == op_frameshift_forward)
			data_.push_back(Packed_operation::frameshift_forward());
		else if (op == op_frameshift_reverse)
			data_.push_back(Packed_operation::frameshift_reverse());
		else if (data_.empty() || data_.back().op() != op || (data_.back().op() == op && data_.back().count() == Packed_operation::MAX_COUNT))
			data_.push_back(Packed_operation(op, 1u));
		else
			++data_.back().code;
	}

	void push_back(Edit_operation op, Letter l)
	{
		data_.push_back(Packed_operation(op, l));
	}

	void push_back(Edit_operation op, unsigned count)
	{
		while (count > 0) {
			const unsigned n = std::min(count, (unsigned)Packed_operation::MAX_COUNT);
			data_.push_back(Packed_operation(op, n));
			count -= n;
		}
	}

	void reverse(size_t begin = 0)
	{
		std::reverse(data_.begin() + begin, data_.end());
	}

	size_t raw_length() const
	{
		return data_.size();
	}

	void push_terminator()
	{
		data_.push_back(Packed_operation::terminator());
	}

	void clear()
	{
		data_.clear();
	}

	void push_back(const Sequence &s, const Reversed&)
	{
		const int l = (int)s.length();
		for (int i = l - 1; i >= 0; --i)
			push_back(op_deletion, s[i]);
	}

	void push_back(const Sequence &s)
	{
		const int l = (int)s.length();
		for (int i = 0; i < l; ++i)
			push_back(op_deletion, s[i]);
	}

	void reserve(size_t n)
	{
		data_.reserve(n);
	}

private:

	std::vector<Packed_operation> data_;

	friend struct Hsp;

};
