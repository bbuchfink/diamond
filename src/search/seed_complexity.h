/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef SEED_COMPLEXITY_H_
#define SEED_COMPLEXITY_H_

#include <math.h>
#include "../basic/value.h"
#include "../basic/shape.h"
#include "../basic/reduction.h"
#include "../basic/config.h"

struct SeedComplexity
{

	static void init(const Reduction &r)
	{
		double p[20];
		for (int i = 0; i < 20; ++i)
			p[i] = 0;
		for (int i = 0; i < 20; ++i)
			p[r(i)] += background_freq[i];
		for (int i = 0; i < 20; ++i)
			prob_[i] = log(p[r(i)]);
		for (int i = 20; i < AMINO_ACID_COUNT; ++i)
			prob_[i] = 1000;
	}

	static bool complex(const char *seq, Shape shape)
	{
		double p = 0;
		for (unsigned i = 0; i < shape.weight_; ++i)
			p += prob_[seq[shape.positions_[i]]];
		return p <= -config.freq_treshold;
	}

private:

	static double prob_[AMINO_ACID_COUNT];

};

#endif