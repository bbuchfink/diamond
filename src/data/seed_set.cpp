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

#include "seed_set.h"

Seed_set::Seed_set(const Sequence_set &seqs, const shape &sh):
	data_((size_t)pow(Reduction::reduction.size(), sh.length_))
{
	uint64_t key;
	const size_t n = seqs.get_length();
	for (size_t i = 0; i < n; ++i) {
		const sequence seq = seqs[i];
		if (seq.length() < sh.length_) continue;
		const size_t j1 = seq.length() - sh.length_ + 1;
		for (size_t j = 0; j < j1; ++j) {
			if (sh.set_seed(key, &seq[j]))
				data_[key] = true;
		}
	}
}