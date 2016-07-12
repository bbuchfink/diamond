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

#ifndef SEQUENCE_SET_H_
#define SEQUENCE_SET_H_

#include <iostream>
#include <string>
#include <algorithm>
#include <queue>
#include "../basic/sequence.h"
#include "string_set.h"
#include "../util/thread.h"
#include "../basic/shape_config.h"

using std::cout;
using std::endl;
using std::pair;

struct Sequence_set : public String_set<'\xff', 1>
{

	Sequence_set()
	{ }

	Sequence_set(Input_stream &file) :
		String_set<'\xff',1>(file)
	{ }

	void print_stats() const
	{
		verbose_stream << "Sequences = " << this->get_length() << ", letters = " << this->letters() << ", average length = " << this->avg_len() << endl;
	}

	pair<size_t, size_t> len_bounds(size_t min_len) const
	{
		const size_t l(this->get_length());
		size_t max = 0, min = std::numeric_limits<size_t>::max();
		for (size_t i = 0; i < l; ++i) {
			max = std::max(this->length(i), max);
			min = this->length(i) >= min_len ? std::min(this->length(i), min) : min;
		}
		return pair<size_t, size_t>(min, max);
	}

	sequence window_infix(size_t offset, unsigned &left) const
	{
		const Letter* begin(this->data(offset));
		unsigned n(0);
		while (*begin != '\xff' && n <= config.window) {
			--begin;
			++n;
		}
		++begin;
		left = config.window + 1 - n;
		const Letter* end(this->data(offset));
		n = 0;
		while (*end != '\xff' && n < config.window) {
			++end;
			++n;
		}
		return sequence(begin, end - begin);
	}

	sequence fixed_window_infix(size_t offset) const
	{
		const Letter* begin(this->data(offset));
		unsigned n(0);
		while (*begin != '\xff' && n <= config.window) {
			--begin;
			++n;
		}
		++begin;
		const Letter* s(this->data(offset - config.window));
		return sequence(s, 2 * config.window, (int)(begin - s));
	}

	vector<size_t> partition(unsigned n_part) const
	{
		vector<size_t> v;
		const size_t l = (this->letters() + n_part - 1) / n_part;
		v.push_back(0);
		for (unsigned i = 0; i < this->get_length();) {
			size_t n = 0;
			while (i < this->get_length() && n < l)
				n += this->length(i++);
			v.push_back(i);
		}
		for (size_t i = v.size(); i < n_part + 1; ++i)
			v.push_back(this->get_length());
		return v;
	}

	size_t reverse_translated_len(size_t i) const
	{
		const size_t j(i - i % 6);
		const size_t l(this->length(j));
		if (this->length(j + 2) == l)
			return l * 3 + 2;
		else if (this->length(j + 1) == l)
			return l * 3 + 1;
		else
			return l * 3;
	}

	size_t avg_len() const
	{
		return this->letters() / this->get_length();
	}

	template <typename _f>
	void enum_seeds(vector<_f> &f) const
	{
		const vector<size_t> p = this->partition((unsigned)f.size());
		Thread_pool threads;
		for (unsigned i = 0; i < f.size(); ++i)
			threads.push_back(launch_thread(enum_seeds_worker<_f>, &f[i], this, (unsigned)p[i], (unsigned)p[i + 1]));
		threads.join_all();
	}

	template<typename _f, typename _entry>
	void enum_seeds_partitioned() const
	{
		vector<Enum_partitioned_callback<_f, _entry> > v;
		v.reserve(config.threads_);
		::partition<unsigned> p(Hashed_seed::p, config.threads_);
		for (unsigned i = 0; i < config.threads_; ++i) {
			Enum_partitioned_callback<_f, _entry> * tmp = new Enum_partitioned_callback<_f, _entry>(p.getMin(i), p.getMax(i));
			v.push_back(*tmp);
			delete tmp;
		}
		enum_seeds(v);
		/*for (unsigned i = 0; i < config.threads_; ++i)
			v[i].flush_queues();*/
	}

	virtual ~Sequence_set()
	{ }

private:

	template<typename _f>
	static void enum_seeds_worker(_f *f, const Sequence_set *seqs, unsigned begin, unsigned end)
	{
		uint64_t key;
		for (unsigned i = begin; i < end; ++i) {
			const sequence seq = (*seqs)[i];
			for (unsigned shape_id = shape_from; shape_id < shape_to; ++shape_id) {
				const shape& sh = shapes[shape_id];
				if (seq.length() < sh.length_) continue;
				for (unsigned j = 0; j < seq.length() - sh.length_ + 1; ++j) {
					if (sh.set_seed(key, &seq[j]))
						(*f)(Hashed_seed(key), seqs->position(i, j), shape_id);
				}
			}
		}
		f->finish();
	}

	template<typename _f, typename _entry>
	struct Enum_partitioned_callback
	{
		Enum_partitioned_callback(unsigned p_begin, unsigned p_end) :
			p_begin(p_begin),
			p_end(p_end)
		{
			memset(counts, 0, sizeof(counts));
		}
		void operator()(Hashed_seed seed, size_t pos, unsigned shape_id)
		{
			const unsigned p = seed.partition();
			buffers[shape_id][p][counts[shape_id][p]++] = _entry(seed, pos);
			if (counts[shape_id][p] == buffer_size)
				flush_buffer(shape_id, p);
			if ((pos & flush_mask) == 0 && shape_id == shape_from)
				flush_queues();
		}
		void flush_buffer(unsigned shape_id, unsigned p)
		{
			mtx[shape_id][p].lock();
			vector<_entry> &q = queues[shape_id][p];
			unsigned &count = counts[shape_id][p];
			const size_t s = q.size();
			q.resize(s + sizeof(_entry)*count);
			memcpy(q.data() + s, buffers[shape_id][p].begin(), sizeof(_entry)*count);
			mtx[shape_id][p].unlock();
			count = 0;
		}
		void finish()
		{
			for (unsigned shape_id = shape_from; shape_id < shape_to;++shape_id)
				for (unsigned p = 0; p < Hashed_seed::p; ++p)
					flush_buffer(shape_id, p);
		}
		void flush_queues()
		{
			for (unsigned shape_id = shape_from; shape_id < shape_to;++shape_id)
				for (unsigned p = p_begin; p < p_end; ++p) {
					mtx[shape_id][p].lock();
					vector<_entry> &q = queues[shape_id][p];
					const size_t size = q.size();
					if (size == 0) {
						mtx[shape_id][p].unlock();
						continue;
					}
					out_buf.resize(size);
					memcpy(out_buf.data(), q.data(), size * sizeof(_entry));
					q.clear();
					mtx[shape_id][p].unlock();
					for (typename vector<_entry>::const_iterator i = out_buf.begin(); i != out_buf.end(); ++i)
						_f()(shape_id, *i);
				}
		}
		enum { buffer_size = 16, flush_mask = 1023 };
		Array<_entry, buffer_size> buffers[Const::max_shapes][Hashed_seed::p];
		unsigned counts[Const::max_shapes][Hashed_seed::p], p_begin, p_end;
		vector<_entry> out_buf;
		static tthread::mutex mtx[Const::max_shapes][Hashed_seed::p];
		static vector<_entry> queues[Const::max_shapes][Hashed_seed::p];
	};

};

template<typename _f, typename _entry>
tthread::mutex Sequence_set::Enum_partitioned_callback<_f, _entry>::mtx[Const::max_shapes][Hashed_seed::p];
template<typename _f, typename _entry>
vector<_entry> Sequence_set::Enum_partitioned_callback<_f, _entry>::queues[Const::max_shapes][Hashed_seed::p];

#endif /* SEQUENCE_SET_H_ */
