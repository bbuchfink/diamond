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

#ifndef PACKED_TRANSCRIPT_H_
#define PACKED_TRANSCRIPT_H_

#include "../util/binary_buffer.h"
#include "../basic/value.h"

typedef enum { op_match=0, op_insertion=1, op_deletion=2, op_substitution=3 } Edit_operation;

struct Packed_operation
{
	Packed_operation(uint8_t code):
		code (code)
	{ }
	Packed_operation(Edit_operation op, unsigned count):
		code ((op<<6) | count)
	{ }
	Packed_operation(Edit_operation op, Letter v):
		code ((op<<6) | (int)v)
	{ }
	operator uint8_t() const
	{ return code; }
	Edit_operation op() const
	{ return (Edit_operation)(code>>6); }
	unsigned count() const
	{
		switch (op()) {
		case op_match:
		case op_insertion:
			return code & 63;
		default:
			return 1;
		}
	}
	Letter letter() const
	{ return code&63; }
	static Packed_operation terminator()
	{ return Packed_operation(op_match, 0u); }
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
			if(op_.op == op_deletion || op_.op == op_substitution) {
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

	void read(Binary_buffer::Iterator &it)
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

	const vector<Packed_operation>& data() const
	{ return data_; }

	const Packed_operation* ptr() const
	{
		return &data_[0];
	}

	void push_back(Edit_operation op)
	{
		if (data_.empty() || data_.back().op() != op || (data_.back().op() == op && data_.back().count() == 63))
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
			const unsigned n = std::min(count, 63u);
			data_.push_back(Packed_operation(op, n));
			count -= n;
		}
	}

	void reverse()
	{
		std::reverse(data_.begin(), data_.end());
	}

	void push_terminator()
	{
		data_.push_back(Packed_operation::terminator());
	}

	void clear()
	{
		data_.clear();
	}

private:

	vector<Packed_operation> data_;

	friend struct Hsp_data;

};

#endif /* PACKED_TRANSCRIPT_H_ */
