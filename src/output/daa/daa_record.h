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
#include <unordered_map>
#include <limits>
#include "daa_file.h"
#include "../basic/packed_sequence.h"
#include "../basic/value.h"
#include "../util/sequence/translate.h"
#include "../basic/packed_transcript.h"
#include "../stats/score_matrix.h"
#include "../basic/match.h"

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

	DAA_query_record(const DAA_file& file, const BinaryBuffer &buf, size_t query_num):
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

	std::string query_name;
	size_t query_num;
	std::vector<Letter> source_seq;
	std::vector<Letter> context[6];
	TranslatedSequence query_seq;

private:

	BinaryBuffer::Iterator init(const BinaryBuffer &buf);

	const DAA_file& file_;
	const BinaryBuffer::Iterator it_;
	
};


void copy_match_record_raw(BinaryBuffer::Iterator& it, TextBuffer& buf, const std::unordered_map<uint32_t, uint32_t>& subject_map);