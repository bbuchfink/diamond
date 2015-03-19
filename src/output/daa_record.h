/****
Copyright (c) 2014, University of Tuebingen
Author: Benjamin Buchfink
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

#include "daa_file.h"
#include "../basic/packed_sequence.h"

using std::string;
using std::vector;

template<typename _val>
struct DAA_match_record;

template<typename _val>
void translate_query(const vector<Nucleotide>& query, vector<_val> *context)
{
	context[0] = query;
	context[1] = Translator::reverse(query);
}

template<>
void translate_query<Amino_acid>(const vector<Nucleotide>& query, vector<Amino_acid> *context)
{
	Translator::translate(query, context);
}

template<typename _val>
struct DAA_query_record
{

	struct Match_iterator
	{
		Match_iterator(const DAA_query_record &parent, Binary_buffer::Iterator it):
			r_ (parent),
			it_ (it),
			good_ (true)
		{ operator++(); }
		const DAA_match_record<_val>& operator*() const
		{ return r_; }
		const DAA_match_record<_val>* operator->() const
		{ return &r_; }
		bool good() const
		{ return good_; }
		Match_iterator& operator++()
		{ if(it_.good()) it_ >> r_; else good_ = false; return *this; }
	private:
		DAA_match_record<_val> r_;
		Binary_buffer::Iterator it_;
		bool good_;
	};

	DAA_query_record(const DAA_file& file, const Binary_buffer &buf):
		file_ (file),
		it_ (init(buf))
	{
		//cout << sequence<Nucleotide>(source_seq.data(), source_seq.size());
	}

	Match_iterator begin() const
	{ return Match_iterator (*this, it_); }

	string query_name;
	vector<Nucleotide> source_seq;
	vector<_val> context[6];

private:

	Binary_buffer::Iterator init(const Binary_buffer &buf)
	{
		Binary_buffer::Iterator it (buf.begin());
		uint32_t query_len;
		it >> query_len;
		it >> query_name;
		uint8_t flags;
		it >> flags;
		if(file_.mode() == blastp) {
			Packed_sequence seq (it, query_len, false, 5);
			seq.unpack(context[0], 5, query_len);
		} else {
			const bool have_n = (flags&1) == 1;
			Packed_sequence seq (it, query_len, have_n, have_n ? 3 : 2);
			seq.unpack(source_seq, have_n ? 3 : 2, query_len);
			translate_query<_val>(source_seq, context);
		}
		return it;
	}

	const DAA_file& file_;
	const Binary_buffer::Iterator it_;

	friend struct DAA_match_record<_val>;
	friend Binary_buffer::Iterator& operator>>(Binary_buffer::Iterator &it, DAA_match_record<_val> &r);

};

template<typename _val>
struct DAA_match_record
{

	DAA_match_record(const DAA_query_record<_val> &query_record):
		parent_ (query_record)
	{ }

	const string& query_name() const
	{ return parent_.query_name; }

	const vector<_val>& query() const
	{ return parent_.context[frame]; }

	size_t db_letters() const
	{ return parent_.file_.db_letters(); }

	unsigned query_end() const
	{
		if(parent_.file_.mode() == blastp) {
			return query_begin + translated_query_len - 1;
		} else if(parent_.file_.mode() == blastx) {
			int len = (int)translated_query_len*3*(frame>2 ? -1 : 1);
			return (int)query_begin + (len > 0 ? -1 : 1) + len;
		} else if(parent_.file_.mode() == blastn) {
			int len = (int)translated_query_len*(frame>0 ? -1 : 1);
			return (int)query_begin + (len > 0 ? -1 : 1) + len;
		} else
			return 0;
	}

	uint32_t total_subject_len, score, query_begin, subject_begin, frame, translated_query_begin, translated_query_len, subject_len, len, identities, mismatches, gap_openings;
	string subject_name;
	Packed_transcript transcript;

private:

	friend Binary_buffer::Iterator& operator>>(Binary_buffer::Iterator &it, DAA_match_record &r)
	{
		uint32_t subject_id;
		it >> subject_id;
		uint8_t flag;
		it >> flag;
		it.read_packed(flag & 3, r.score);
		it.read_packed((flag>>2)&3, r.query_begin);
		it.read_packed((flag>>4)&3, r.subject_begin);
		r.transcript.read(it);
		r.subject_name = r.parent_.file_.ref_name(subject_id);
		r.total_subject_len = r.parent_.file_.ref_len(subject_id);
		if(r.parent_.file_.mode() == blastx) {
			r.frame = (flag&(1<<6)) == 0 ? r.query_begin % 3 : 3+(r.parent_.source_seq.size() - 1 - r.query_begin)%3;
			r.translated_query_begin = query_translated_begin<_val>(r.query_begin, r.frame, r.parent_.source_seq.size(), true);
		} else if (r.parent_.file_.mode() == blastp) {
			r.frame = 0;
			r.translated_query_begin = r.query_begin;
		} else {
			r.frame = (flag&(1<<6)) == 0 ? 0 : 1;
			r.translated_query_begin = query_translated_begin<_val>(r.query_begin, r.frame, r.parent_.source_seq.size(), false);
		}
		r.parse();
		return it;
	}

	void parse()
	{
		translated_query_len = 0;
		subject_len = 0;
		len = 0;
		identities = 0;
		mismatches = 0;
		gap_openings = 0;
		unsigned d = 0;
		for(Packed_transcript::Const_iterator<char> i = transcript.template begin<char>(); i.good(); ++i) {
			len += i->count;
			switch(i->op) {
			case op_match:
				identities += i->count;
				translated_query_len += i->count;
				subject_len += i->count;
				d = 0;
				break;
			case op_substitution:
				mismatches += i->count;
				translated_query_len += i->count;
				subject_len += i->count;
				d = 0;
				break;
			case op_insertion:
				translated_query_len += i->count;
				++gap_openings;
				d = 0;
				break;
			case op_deletion:
				subject_len += i->count;
				if(d == 0)
					++gap_openings;
				d += i->count;
			}
		}
	}

	const DAA_query_record<_val> &parent_;

};

#endif /* DAA_RECORD_H_ */
