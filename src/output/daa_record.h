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

	struct Match
	{

		Match(const DAA_query_record &query_record) :
			hit_num(std::numeric_limits<uint32_t>::max()),
			subject_id(std::numeric_limits<uint32_t>::max()),
			parent_(query_record)
		{ }

		const string& query_name() const;
		const vector<Letter>& query() const;
		uint64_t db_letters() const;
		unsigned query_end() const;

		double bit_score() const
		{
			return score_matrix.bitscore(score);
		}

		double evalue() const
		{
			return score_matrix.evalue(score, (size_t)db_letters(), (unsigned)query().size());
		}

		bool is_query_translated() const
		{
			return parent_.file_.mode() == mode_blastx;
		}

		int blast_query_frame() const
		{
			return is_query_translated() ? (frame <= 2 ? (int)frame + 1 : 2 - (int)frame) : 0;
		}

		std::pair<unsigned, unsigned> unoriented_query_range() const
		{
			if (frame > 2)
				return std::pair<unsigned, unsigned>(query_end(), query_begin);
			else
				return std::pair<unsigned, unsigned>(query_begin, query_end());
		}

		struct Position_iterator
		{
			Position_iterator(const Match &parent) :
				query_pos(parent.translated_query_begin),
				subject_pos(parent.subject_begin),
				ptr_(parent.transcript.ptr()),
				count_(ptr_->count()),
				parent_(parent)
			{ }
			bool good() const
			{
				return *ptr_ != Packed_operation::terminator();
			}
			Position_iterator& operator++()
			{
				switch (op()) {
				case op_deletion:
					++subject_pos;
					break;
				case op_insertion:
					++query_pos;
					break;
				case op_match:
				case op_substitution:
					++query_pos;
					++subject_pos;
				}
				--count_;
				if (count_ == 0) {
					++ptr_;
					count_ = ptr_->count();
				}
				return *this;
			}
			Edit_operation op() const
			{
				return ptr_->op();
			}
			Letter query() const
			{
				return parent_.query()[query_pos];
			}
			Letter subject() const
			{
				switch (op()) {
				case op_substitution:
				case op_deletion:
					return ptr_->letter();
				default:
					return query();
				}				
			}
			char query_char() const
			{
				switch (op()) {
				case op_deletion:
					return '-';
				default:
					return value_traits.alphabet[(long)query()];
				}
			}
			char subject_char() const
			{
				switch (op()) {
				case op_insertion:
					return '-';
				default:
					return value_traits.alphabet[(long)subject()];
				}
			}
			char midline_char() const
			{
				switch (op()) {
				case op_match:
					return value_traits.alphabet[(long)query()];
				case op_substitution:
					return score() > 0 ? '+' : ' ';
				default:
					return ' ';
				}
			}
			int score() const
			{
				return score_matrix(query(), subject());
			}
			unsigned query_pos, subject_pos;
		private:
			const Packed_operation *ptr_;
			unsigned count_;
			const Match &parent_;
		};
		
		Position_iterator begin() const
		{
			return Position_iterator(*this);
		}

		uint32_t total_subject_len, score, query_begin, subject_begin, frame, translated_query_begin, translated_query_len, subject_len, len, identities, mismatches, positives, gap_openings, gaps, hsp_num, hit_num, subject_id;
		string subject_name;
		Packed_transcript transcript;

	private:

		void parse();
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
		const Match& operator*() const
		{
			return r_;
		}
		const Match* operator->() const
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
		switch (file_.mode()) {
		case mode_blastx:
			return source_seq.size();
		default:
			return context[0].size();
		}
	}

	size_t db_letters() const
	{
		return file_.db_letters();
	}

	size_t db_seqs() const
	{
		return file_.db_seqs();
	}

	double lambda() const
	{
		return file_.lambda();
	}

	double kappa() const
	{
		return file_.kappa();
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
