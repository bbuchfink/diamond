/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
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

#ifndef DAA_RECORD_H_
#define DAA_RECORD_H_

#include <limits>
#include "daa_file.h"
#include "../basic/packed_sequence.h"
#include "../basic/value.h"
#include "../basic/translate.h"
#include "../basic/packed_transcript.h"
#include "../align/match_func.h"
#include "../basic/score_matrix.h"

using std::string;
using std::vector;

/*void translate_query(const vector<Letter>& query, vector<Letter> *context)
{
	context[0] = query;
	context[1] = Translator::reverse(query);
}*/

inline void translate_query(const vector<Letter>& query, vector<Letter> *context)
{
	Translator::translate(query, context);
}

struct DAA_query_record
{

	struct Match : public Hsp_data
	{

		Match(const DAA_query_record &query_record) :
			hit_num(std::numeric_limits<uint32_t>::max()),
			subject_id(std::numeric_limits<uint32_t>::max()),
			parent_(query_record)
		{ }

		Hsp_context context()
		{
			return Hsp_context(*this,
				(unsigned)parent_.query_num,
				sequence(parent_.context[frame]),
				align_mode.query_translated ? sequence(parent_.source_seq) : sequence(parent_.context[0]),
				parent_.query_name.c_str(),
				subject_id,
				subject_id,
				subject_name.c_str(),
				subject_len,
				hit_num,
				hsp_num);
		}

		uint32_t hsp_num, hit_num, subject_id, subject_len;
		string subject_name;

	private:

		const DAA_query_record &parent_;
		friend Binary_buffer::Iterator& operator>>(Binary_buffer::Iterator &it, Match &r);

	};

	struct Match_iterator
	{
		Match_iterator(const DAA_query_record &parent, Binary_buffer::Iterator it) :
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
			if (it_.good()) it_ >> r_; else good_ = false; return *this;
		}
	private:
		Match r_;
		Binary_buffer::Iterator it_;
		bool good_;
	};

	DAA_query_record(const DAA_file& file, const Binary_buffer &buf, size_t query_num):
		query_num (query_num),
		file_(file),
		it_(init(buf))
	{ }

	Match_iterator begin() const
	{
		return Match_iterator(*this, it_);
	}

	size_t query_len() const
	{
		return align_mode.query_translated ? source_seq.size() : context[0].size();
	}

	string query_name;
	size_t query_num;
	vector<Letter> source_seq;
	vector<Letter> context[6];

private:

	Binary_buffer::Iterator init(const Binary_buffer &buf);

	const DAA_file& file_;
	const Binary_buffer::Iterator it_;
	
	friend Binary_buffer::Iterator& operator>>(Binary_buffer::Iterator &it, Match &r);

};

Binary_buffer::Iterator& operator>>(Binary_buffer::Iterator &it, DAA_query_record::Match &r);

#endif /* DAA_RECORD_H_ */
