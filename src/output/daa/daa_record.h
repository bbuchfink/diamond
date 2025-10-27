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
#include <unordered_map>
#include <limits>
#include "daa_file.h"
#include "basic/value.h"
#include "util/sequence/translate.h"
#include "basic/match.h"

inline void translate_query(const std::vector<Letter>& query, std::vector<Letter> *context)
{
	Translator::translate(query, context);
}

struct DAA_query_record
{

	struct Match : public Hsp
	{

		Match(const DAA_query_record &query_record) :
			Hsp(true),
			hit_num(std::numeric_limits<uint32_t>::max()),
			subject_id(std::numeric_limits<uint32_t>::max()),
			parent_(query_record)
		{ }

		HspContext context()
		{
			return HspContext(*this,
				(unsigned)parent_.query_num,
				0,
				parent_.query_seq,
				parent_.query_name.c_str(),
				subject_id,
				subject_len,
				subject_name.c_str(),
				hit_num,
				hsp_num,
				Sequence());
		}

		void read(BinaryBuffer::Iterator &it);

		uint32_t hsp_num, hit_num, subject_id, subject_len;
		std::string subject_name;

	private:

		const DAA_query_record &parent_;

	};

	struct Match_iterator
	{
		Match_iterator(const DAA_query_record &parent, BinaryBuffer::Iterator it) :
			r_(parent),
			it_(it),
			good_(true)
		{
			operator++();
		}
		Match& operator*()
		{
			return r_;
		}
		Match* operator->()
		{
			return &r_;
		}
		bool good() const
		{
			return good_;
		}
		Match_iterator& operator++()
		{
			if (it_.good()) r_.read(it_); else good_ = false; return *this;
		}
	private:
		Match r_;
		BinaryBuffer::Iterator it_;
		bool good_;
	};

	DAA_query_record(const DAAFile& file, const BinaryBuffer &buf, size_t query_num):
		query_num (query_num),
		file_(file),
		it_(init(buf))
	{ }

	Match_iterator begin() const
	{
		return Match_iterator(*this, it_);
	}

	BinaryBuffer::Iterator raw_begin() const {
		return it_;
	}

	size_t query_len() const
	{
		return align_mode.query_translated ? source_seq.size() : context[0].size();
	}

	const DAAFile& file() const {
		return file_;
	}

	std::string query_name;
	size_t query_num;
	std::vector<Letter> source_seq;
	std::vector<Letter> context[6];
	TranslatedSequence query_seq;

private:

	BinaryBuffer::Iterator init(const BinaryBuffer &buf);

	const DAAFile& file_;
	const BinaryBuffer::Iterator it_;
	
};


void copy_match_record_raw(BinaryBuffer::Iterator& it, TextBuffer& buf, const std::unordered_map<uint32_t, uint32_t>& subject_map);