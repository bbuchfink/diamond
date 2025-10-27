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

#include <atomic>
#include <thread>
#include "sequence_set.h"
#include "util/sequence/sequence.h"
#include "util/log_stream.h"

using std::vector;
using std::pair;
using std::string;

SequenceSet::SequenceSet(Alphabet alphabet) :
	alphabet_(alphabet)
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
	return std::pair<Length, Length>(min, max);
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

void SequenceSet::convert_to_std_alph(size_t id)
{
	if (alphabet_ == Alphabet::STD)
		return;
	Letter* ptr = this->ptr(id);
	const size_t len = length(id);
	alph_ncbi_to_std(ptr, ptr + len);
}

void SequenceSet::convert_all_to_std_alph(size_t threads)
{
	if (alphabet_ == Alphabet::STD)
		return;
	std::atomic_size_t next(0);
	auto worker = [this, &next] {
		const size_t n = this->size();
		size_t i;
		while ((i = next++) < n)
			this->convert_to_std_alph(i);
	};
	vector<std::thread> t;
	for (size_t i = 0; i < threads; ++i)
		t.emplace_back(worker);
	for (auto& i : t)
		i.join();
	alphabet_ = Alphabet::STD;
}

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