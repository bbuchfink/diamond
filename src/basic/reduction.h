/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include <vector>
#include <array>
#include <string>
#include "value.h"
#include "sequence.h"
#include "util/memory_resource.h"

struct Reduction
{

	Reduction(const char* definition_string);

	unsigned size() const
	{
		return size_;
	}

	int bit_size() const
	{
		return bit_size_;
	}

	double bit_size_exact() const {
		return bit_size_exact_;
	}

	unsigned operator()(Letter a) const
	{
		return map_[(long)a];
	}

	unsigned operator()(size_t a) const
	{
		return map_[a];
	}

	const Letter* map8() const
	{
		return map8_;
	}

	const Letter* map8b() const
	{
		return map8b_;
	}

	inline friend std::ostream& operator<<(std::ostream &os, const Reduction &r)
	{
		for (unsigned i = 0; i < r.size_; ++i) {
			os << '[';
			for (unsigned j = 0; j < 20; ++j)
				if (r.map_[j] == i)
					os << value_traits.alphabet[j];
			os << ']';
		}
		return os;
	}

	static Reduction& get_reduction()
	{
		return instance;
	}

	static void set_reduction(const Reduction& r)
	{
		get_reduction() = r;
	}

	static std::pmr::vector<Letter> reduce_seq(const Sequence &seq, std::pmr::monotonic_buffer_resource& pool)
	{
		std::pmr::vector<Letter> dst(&pool);
		dst.reserve(seq.length());
		const Reduction& reduction_instance = get_reduction();
		for (Loc i = 0; i < seq.length(); ++i)
			dst.push_back(reduction_instance(seq[i]));
		return dst;
	}

	double freq(unsigned bucket) const {
		return freq_[bucket];
	}

	std::string decode_seed(const uint64_t seed, const size_t len) const;

private:

	unsigned map_[256];
	alignas(16) Letter map8_[256], map8b_[256];
	unsigned size_;
	int bit_size_;
	double bit_size_exact_;
	std::array<double, TRUE_AA> freq_;
	static Reduction instance;

};