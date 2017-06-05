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

#ifndef TRACE_PT_BUFFER_H_
#define TRACE_PT_BUFFER_H_

#include "../util/async_buffer.h"
#include "../basic/match.h"

#pragma pack(1)

struct hit
{
	typedef uint32_t Seed_offset;

	unsigned	query_;
	Packed_loc	subject_;
	Seed_offset	seed_offset_;
	hit() :
		query_(),
		subject_(),
		seed_offset_()
	{ }
	hit(unsigned query, Packed_loc subject, Seed_offset seed_offset) :
		query_(query),
		subject_(subject),
		seed_offset_(seed_offset)
	{ }
	bool operator<(const hit &rhs) const
	{
		return query_ < rhs.query_;
	}
	bool blank() const
	{
		return subject_ == 0;
	}
	unsigned operator%(unsigned i) const
	{
		return (query_ / align_mode.query_contexts) % i;
	}
	unsigned operator/(size_t i) const
	{
		return (query_ / align_mode.query_contexts) / (unsigned)i;
	}
	int64_t global_diagonal() const
	{
		return (int64_t)subject_ - (int64_t)seed_offset_;
	}
	template<unsigned _d>
	static unsigned query_id(const hit& x)
	{
		return x.query_ / _d;
	}
	template<unsigned _d>
	struct Query_id
	{
		unsigned operator()(const hit& x) const
		{
			return query_id<_d>(x);
		}
	};
	static bool cmp_subject(const hit &lhs, const hit &rhs)
	{
		return lhs.subject_ < rhs.subject_
			|| (lhs.subject_ == rhs.subject_ && lhs.seed_offset_ < rhs.seed_offset_);
	}
	static bool cmp_normalized_subject(const hit &lhs, const hit &rhs)
	{
		const uint64_t x = (uint64_t)lhs.subject_ + (uint64_t)rhs.seed_offset_, y = (uint64_t)rhs.subject_ + (uint64_t)lhs.seed_offset_;
		return x < y || (x == y && lhs.seed_offset_ < rhs.seed_offset_);
	}
	friend std::ostream& operator<<(std::ostream &s, const hit &me)
	{
		s << me.query_ << '\t' << me.subject_ << '\t' << me.seed_offset_ << '\n';
		return s;
	}
} PACKED_ATTRIBUTE;

#pragma pack()

struct Trace_pt_buffer : public Async_buffer<hit>
{
	Trace_pt_buffer(size_t input_size, const string &tmpdir, unsigned query_bins):
		Async_buffer<hit>(input_size, tmpdir, query_bins)
	{}
	static Trace_pt_buffer *instance;
};

struct Trace_pt_list : public vector<hit>
{
	void init()
	{
		pos_ = this->begin();
		total_ = 0;
		count_ = 1;
#ifdef PRE_PARTITION
		p_.clear();
		p_.push_back(0);
		idx_ = 0;
		const unsigned c = query_contexts();
		typename vector<hit<_locr,_locl> >::iterator i = this->begin();
		unsigned total=0,count=1;
		for(; i < this->end();) {
			unsigned n=0;
			const unsigned min_size = std::max(4*total/count/5 + 1, program_options::fetch_size);
			for(;i<this->end() && n<min_size;) {
				const unsigned q = i->query_/c;
				for(; i<this->end() && i->query_/c == q; ++i)
					++n;
			}
			++count;
			total += n;
			p_.push_back(i - this->begin());
		}
		p_.push_back(i - this->begin());
#endif
	}
	struct Query_range
	{
		Query_range(Trace_pt_list &parent):
			parent_ (parent)
		{ }
#ifndef PRE_PARTITION
		bool operator()()
		{

			begin = parent_.pos_;
			//end = std::min(std::max(begin + 3*parent_.total_/parent_.count_/4 + 1, begin+program_options::fetch_size), parent_.end());
#ifdef NDEBUG
			end = std::min(begin + 3*parent_.total_/parent_.count_/4 + 1, parent_.end());
#else
			//end = parent_.end();
			ptrdiff_t x = std::min((ptrdiff_t)(3 * parent_.total_ / parent_.count_ / 4 + 1), parent_.end() - begin);
			end = std::min(begin + x, parent_.end());
#endif
			if(end >= parent_.end())
				return false;
			const unsigned c = align_mode.query_contexts, q = end->query_/c;
			for(; end<parent_.end() && end->query_/c == q; ++end);
			parent_.pos_ = end;
			parent_.total_ += end - begin;
			++parent_.count_;
			return end < parent_.end();
		}
#else
		bool operator()()
		{
			begin = parent_.begin()+parent_.p_[parent_.idx_];
			end = parent_.begin()+parent_.p_[parent_.idx_+1];
			printf("%lu %lu %lu\n", parent_.p_[parent_.idx_], parent_.p_[parent_.idx_+1], parent_.p_[parent_.idx_+1]-parent_.p_[parent_.idx_]);
			++parent_.idx_;
			return parent_.idx_ < parent_.p_.size()-1;
		}
#endif
		Trace_pt_list::iterator begin, end;
	private:
		Trace_pt_list &parent_;
	};
	Query_range get_range()
	{ return Query_range (*this); }
private:
	vector<hit>::iterator pos_;
#ifdef PRE_PARTITION
	vector<size_t> p_;
	unsigned idx_;
#else
	size_t total_, count_;
#endif
};

#endif /* TRACE_PT_BUFFER_H_ */

