/****
Copyright (c) 2015-2016, University of Tuebingen, Benjamin Buchfink
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

	void read(Buffered_file &f)
	{
		data_.clear();
		uint8_t code;
		do {
			f.read(code);
			data_.push_back(code);
		} while (code != Packed_operation::terminator());
	}

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
		data_.push_back(Packed_operation(op, count));
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
