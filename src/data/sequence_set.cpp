/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include <atomic>
#include <thread>
#include "sequence_set.h"
#include "util/sequence/sequence.h"
#include "util/log_stream.h"

using std::vector;
using std::pair;
using std::string;

SequenceSet::SequenceSet()
{ }

void SequenceSet::print_stats() const
{
	verbose_stream << "Sequences = " << this->size() << ", letters = " << this->letters() << ", average length = " << this->avg_len() << std::endl;
}

std::pair<SequenceSet::Length, SequenceSet::Length> SequenceSet::len_bounds(Length min_len) const
{
	const size_t l(this->size());
	Length max = 0, min = std::numeric_limits<Length>::max();
	for (size_t i = 0; i < l; ++i) {
		max = std::max(this->length(i), max);
		min = this->length(i) >= min_len ? std::min(this->length(i), min) : min;
	}
	return pair<Length, Length>(min, max);
}

SequenceSet::Length SequenceSet::max_len(size_t begin, size_t end) const
{
	Length max = 0;
	for (size_t i = begin; i < end; ++i)
		max = std::max(this->length(i), max);
	return max;
}

std::vector<uint32_t> SequenceSet::partition(unsigned n_part, bool shortened, bool context_reduced) const
{
	std::vector<uint32_t> v;
	const size_t l = (this->letters() + n_part - 1) / n_part;
	const Id c = context_reduced ? align_mode.query_contexts : 1;
	if(!shortened)
		v.push_back(0);
	for (Id i = 0; i < this->size();) {
		size_t n = 0;
		while (i < this->size() && n < l) {
			for (Id j = 0; j < c; ++j)
				n += this->length(i++);
		}
		v.push_back(i / c);
	}
	for (size_t i = v.size(); i < n_part + (shortened ? 0 : 1); ++i)
		v.push_back(this->size() / c);
	return v;
}

size_t SequenceSet::reverse_translated_len(size_t i) const
{
	const size_t j(i - i % 6);
	const Loc l(this->length(j));
	if (this->length(j + 2) == l)
		return l * 3 + 2;
	else if (this->length(j + 1) == l)
		return l * 3 + 1;
	else
		return l * 3;
}

TranslatedSequence SequenceSet::translated_seq(const Sequence& source, size_t i) const
{
	if (!align_mode.query_translated)
		return TranslatedSequence((*this)[i]);
	return TranslatedSequence(source, (*this)[i], (*this)[i + 1], (*this)[i + 2], (*this)[i + 3], (*this)[i + 4], (*this)[i + 5]);
}

size_t SequenceSet::avg_len() const
{
	return this->letters() / this->size();
}

SequenceSet::~SequenceSet()
{ }

size_t max_id_len(const StringSet& ids)
{
	size_t max(0);
	for (BlockId i = 0; i < ids.size(); ++i)
		max = std::max(max, find_first_of(ids[i], Util::Seq::id_delimiters));
	return max;
}

std::vector<std::pair<Loc, BlockId>> SequenceSet::lengths() const {
	vector<pair<Loc, BlockId>> l;
	l.reserve(size());
	for (BlockId i = 0; i < size(); ++i)
		l.emplace_back(length(i), i);
	return l;
}

Loc SequenceSet::source_length(BlockId i) const {
	if (align_mode.query_contexts == 1)
		return length(i);
	const BlockId j = i - i % align_mode.query_contexts;
	return length(j) + length(j + 1) + length(j + 2) + 2;
}