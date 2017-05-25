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
	void enum_seeds(vector<_f> &f, const vector<size_t> &p, size_t shape_begin, size_t shape_end) const
	{
		Thread_pool threads;
		for (unsigned i = 0; i < f.size(); ++i)
			threads.push_back(launch_thread(enum_seeds_worker<_f>, &f[i], this, (unsigned)p[i], (unsigned)p[i + 1]));
		threads.join_all();
	}

	virtual ~Sequence_set()
	{ }

private:

	template<typename _f>
	static void enum_seeds_worker(_f *f, const Sequence_set *seqs, unsigned begin, unsigned end, size_t shape_begin, size_t shape_end)
	{
		uint64_t key;
		for (unsigned i = begin; i < end; ++i) {
			const sequence seq = (*seqs)[i];
			for (unsigned shape_id = shape_begin; shape_id < shape_end ++shape_id) {
				const shape& sh = shapes[shape_id];
				if (seq.length() < sh.length_) continue;
				for (unsigned j = 0; j < seq.length() - sh.length_ + 1; ++j) {
					if (sh.set_seed(key, &seq[j]))
						(*f)(key, seqs->position(i, j), shape_id);
				}
			}
		}
		f->finish();
	}

};

#endif /* SEQUENCE_SET_H_ */