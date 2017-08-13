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

#ifndef MATCH_H_
#define MATCH_H_

#include <limits>
#include "sequence.h"
#include "../util/util.h"
#include "../util/async_buffer.h"
#include "packed_loc.h"
#include "../util/system.h"
#include "value.h"
#include "packed_transcript.h"
#include "score_matrix.h"

enum Strand { FORWARD, REVERSE };

inline interval normalized_range(unsigned pos, int len, Strand strand)
{
	return strand == FORWARD
			? interval (pos, pos + len)
			: interval (pos + 1 + len, pos + 1);
}

interval untranslate_range(const interval &r, unsigned frame, size_t l);

struct Diagonal_segment
{
	Diagonal_segment():
		len(0)
	{}
	Diagonal_segment(int query_pos, int subject_pos, int len, int score):
		i(query_pos),
		j(subject_pos),
		len(len),
		score (score)
	{}
	bool empty() const
	{
		return len == 0;
	}
	interval query_range() const
	{
		return interval(i, i + len);
	}
	interval subject_range() const
	{
		return interval(j, j + len);
	}
	int subject_last() const
	{
		return j + len - 1;
	}
	int query_last() const
	{
		return i + len - 1;
	}
	int subject_end() const
	{
		return j + len;
	}
	int query_end() const
	{
		return i + len;
	}
	int diag() const
	{
		return i - j;
	}
	Diagonal_segment intersect(const Diagonal_segment &x) const
	{
		if (diag() != x.diag())
			return Diagonal_segment();
		else {
			const interval q = ::intersect(query_range(), x.query_range());
			return Diagonal_segment(q.begin_, ::intersect(subject_range(), x.subject_range()).begin_, q.length(), 0);
		}
	}
	bool is_enveloped(const Diagonal_segment &x) const
	{
		return score <= x.score
			&& query_range().overlap_factor(x.query_range()) == 1
			&& subject_range().overlap_factor(x.subject_range()) == 1;
	}
	Diagonal_segment transpose() const
	{
		return Diagonal_segment(j, i, len, score);
	}
	int partial_score(int diff) const
	{
		return score*std::max(len - diff, 0) / len;
	}
	bool operator<=(const Diagonal_segment &rhs) const
	{
		return i + len <= rhs.i && j + len <= rhs.j;
	}
	bool operator==(const Diagonal_segment &rhs) const
	{
		return i == rhs.i && j == rhs.j && len == rhs.len;
	}
	static bool cmp_subject(const Diagonal_segment &x, const Diagonal_segment &y)
	{
		return x.j < y.j || (x.j == y.j && x.i < y.i);
	}
	static bool cmp_subject_end(const Diagonal_segment &x, const Diagonal_segment &y)
	{
		return x.subject_end() < y.subject_end();
	}
	static bool cmp_heuristic(const Diagonal_segment &x, const Diagonal_segment &y)
	{
		return (x.subject_end() < y.subject_end() && x.j < y.j)
			|| (x.j - y.j < y.subject_end() - x.subject_end());
	}
	friend int abs_shift(const Diagonal_segment &x, const Diagonal_segment &y)
	{
		return abs(x.diag() - y.diag());
	}
	friend std::ostream& operator<<(std::ostream &s, const Diagonal_segment &d)
	{
		s << "i=" << d.i << " j=" << d.j << " l=" << d.len << " score=" << d.score;
		return s;
	}
	int i, j, len, score;
};

struct Intermediate_record;

struct Hsp_data
{
	Hsp_data() :
		score(0),
		frame(0),
		length(0),
		identities(0),
		mismatches(0),
		positives(0),
		gap_openings(0),
		gaps(0),
		sw_score(0)
	{}
	Hsp_data(int score) :
		score(unsigned(score)),
		frame(0),
		length(0),
		identities(0),
		mismatches(0),
		positives(0),
		gap_openings(0),
		gaps(0),
		sw_score(0)
	{}
	Hsp_data(const Intermediate_record &r, unsigned query_source_len);
	struct Iterator
	{
		Iterator(const Hsp_data &parent):
			query_pos(parent.query_range.begin_),
			subject_pos(parent.subject_range.begin_),
			ptr_(parent.transcript.ptr()),
			count_(ptr_->count())
		{ }
		bool good() const
		{
			return *ptr_ != Packed_operation::terminator();
		}
		Iterator& operator++()
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
		bool to_next_match()
		{
			do {
				this->operator++();
			} while (good() && (op() == op_deletion || op() == op_insertion));
			return good();
		}
		unsigned query_pos, subject_pos;
	protected:
		const Packed_operation *ptr_;
		unsigned count_;
	};
	Iterator begin() const
	{
		return Iterator(*this);
	}
	void set_source_range(unsigned frame, unsigned dna_len);
	void set_source_range(unsigned oriented_query_begin)
	{
		query_source_range = frame < 3 ? interval(oriented_query_begin, oriented_query_begin + 3 * query_range.length()) : interval(oriented_query_begin + 1 - 3 * query_range.length(), oriented_query_begin + 1);
	}
	interval oriented_range() const
	{
		if (frame < 3)
			return interval(query_source_range.begin_, query_source_range.end_ - 1);
		else
			return interval(query_source_range.end_ - 1, query_source_range.begin_);
	}
	void set_translated_query_begin(unsigned oriented_query_begin, unsigned dna_len)
	{
		int f = frame <= 2 ? frame + 1 : 2 - frame;
		if (f > 0)
			query_range.begin_ = (oriented_query_begin - (f - 1)) / 3;
		else
			query_range.begin_ = (dna_len + f - oriented_query_begin) / 3;
	}
	int blast_query_frame() const
	{
		return align_mode.query_translated ? (frame <= 2 ? (int)frame + 1 : 2 - (int)frame) : 0;
	}
	bool operator<(const Hsp_data &rhs) const
	{
		return score > rhs.score;
	}
	double id_percent() const
	{
		return (double)identities * 100.0 / (double)length;
	}
	double query_cover_percent(unsigned query_source_len) const
	{
		return (double)query_source_range.length() * 100 / query_source_len;
	}
	double subject_cover_percent(unsigned subject_len) const
	{
		return (double)subject_range.length() * 100 / subject_len;
	}
	static bool cmp_query_pos(const Hsp_data &x, const Hsp_data &y)
	{
		return x.query_range.begin_ < y.query_range.begin_;
	}
	bool pass_through(const Diagonal_segment &d) const;
	bool is_weakly_enveloped(const Hsp_data &j) const;
	std::pair<int, int> diagonal_bounds() const;
	void merge(const Hsp_data &right, const Hsp_data &left, unsigned query_anchor, unsigned subject_anchor);
	unsigned score, frame, length, identities, mismatches, positives, gap_openings, gaps, sw_score;
	float time;
	interval query_source_range, query_range, subject_range;
	Packed_transcript transcript;
};

struct Hsp_context
{
	Hsp_context(Hsp_data& hsp, unsigned query_id, const sequence &query, const sequence &source_query, const char *query_name, unsigned subject_id, unsigned orig_subject_id, const char *subject_name, unsigned subject_len, unsigned hit_num, unsigned hsp_num) :		
		query(query),
		source_query(source_query),
		query_name(query_name),
		subject_name(subject_name),
		query_id(query_id),
		subject_id(subject_id),
		orig_subject_id(orig_subject_id),
		subject_len(subject_len),
		hit_num(hit_num),
		hsp_num(hsp_num),
		hsp_(hsp)
	{}
	struct Iterator : public Hsp_data::Iterator
	{
		Iterator(const Hsp_context &parent) :
			Hsp_data::Iterator(parent.hsp_),
			parent_(parent)
		{ }
		Letter query() const
		{
			return parent_.query[query_pos];
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
	private:
		const Hsp_context &parent_;
	};
	Iterator begin() const
	{
		return Iterator(*this);
	}
	Packed_transcript::Const_iterator begin_old() const
	{
		return hsp_.transcript.begin();
	}
	unsigned score() const
	{ return hsp_.score; }
	double evalue() const
	{
		return score_matrix.evalue(score(), config.db_size, (unsigned)query.length());
	}
	double bit_score() const
	{
		return score_matrix.bitscore(score());
	}
	double sw_score() const
	{
		return score_matrix.bitscore(hsp_.sw_score);
	}
	double time() const
	{
		return hsp_.time;
	}
	unsigned frame() const
	{ return hsp_.frame; }
	unsigned length() const
	{ return hsp_.length; }
	unsigned identities() const
	{ return hsp_.identities; }
	unsigned mismatches() const
	{ return hsp_.mismatches; }
	unsigned positives() const
	{ return hsp_.positives; }
	unsigned gap_openings() const
	{ return hsp_.gap_openings; }
	unsigned gaps() const
	{ return hsp_.gaps; }
	const interval& query_source_range() const
	{ return hsp_.query_source_range; }
	const interval& query_range() const
	{ return hsp_.query_range; }
	const interval& subject_range() const
	{ return hsp_.subject_range; }
	interval oriented_query_range() const
	{ return hsp_.oriented_range(); }
	int blast_query_frame() const
	{ return hsp_.blast_query_frame(); }
	Packed_transcript transcript() const
	{ return hsp_.transcript; }
	Hsp_context& parse();
	Hsp_context& set_query_source_range(unsigned oriented_query_begin);

	const sequence query, source_query;
	const char *query_name, *subject_name;
	const unsigned query_id, subject_id, orig_subject_id, subject_len, hit_num, hsp_num;
private:	
	Hsp_data &hsp_;
};

#endif /* MATCH_H_ */
