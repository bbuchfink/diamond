/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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

/* Based on tantan by Martin C. Frith. See:
http://cbrc3.cbrc.jp/~martin/tantan/
A new repeat-masking method enables specific detection of homologous sequences, MC Frith, Nucleic Acids Research 2011 39(4):e23. */

#include <array>
#include <stdint.h>
#include <algorithm>
#include <Eigen/Core>
#include "../basic/value.h"
#include "def.h"

using Eigen::Array;
using Eigen::Dynamic;

namespace Util { namespace tantan { namespace DISPATCH_ARCH {

Mask::Ranges mask(Letter *seq,
	int len,
	const float **likelihood_ratio_matrix,
	float p_repeat,
	float p_repeat_end,
	float repeat_growth,
	float p_mask,
	const Letter *mask_table) {
	constexpr int WINDOW = 50, RESERVE = 50000;

	if (len == 0)
		return {};

	thread_local std::array<Array<float, Dynamic, 1>, AMINO_ACID_COUNT> e;
	thread_local Array<float, Dynamic, 1> pb;
	thread_local Array<float, Dynamic, 1> scale;
	Array<float, WINDOW, 1> f(0.0), d, t;
	const float b2b = 1.0f - p_repeat, f2f = 1.0f - p_repeat_end, b2f0 = p_repeat * (1.0f - repeat_growth) / (1.0f - pow(repeat_growth, (float)WINDOW));
	float b = 1.0;
	
	d[WINDOW - 1] = b2f0;
	for (int i = WINDOW - 2; i >= 0; --i)
		d[i] = d[i + 1] * repeat_growth;

	pb.resize(std::max(len, RESERVE));
	scale.resize(std::max((len - 1) / 16 + 1, (RESERVE - 1) / 16 + 1));

	for (int i = 0; i < (int)AMINO_ACID_COUNT; ++i) {
		e[i].resize(std::max(RESERVE, len + WINDOW));
		const float *l = likelihood_ratio_matrix[i];
		float* p = &e[i][len - 1];
		for (int j = 0; j < len; ++j)
			*(p--) = l[(size_t)letter_mask(seq[j])];
		std::fill(e[i].data() + len, e[i].data() + len + WINDOW, (float)0.0);
	}	

	for (int i = 0; i < len; ++i) {
		const float s = f.sum();
		t = b;
		t *= d;
		f *= f2f;
		f += t;
		
		f *= e[(size_t)letter_mask(seq[i])].template segment<WINDOW>(len - i, WINDOW);
		b = b * b2b + s * p_repeat_end;

		if ((i & 15) == 15) {
			const float s = 1 / b;
			scale[i / 16] = s;
			b *= s;
			f *= s;
		}

		pb[i] = b;
	}

	const float z = b * b2b + f.sum() * p_repeat_end;
	b = b2b;
	f = p_repeat_end;
	Mask::Ranges ranges;

	for (int i = len - 1; i >= 0; --i) {
		const float pf = 1 - (pb[i] * b / z);

		if ((i & 15) == 15) {
			const float s = scale[i / 16];
			b *= s;
			f *= s;
		}

		f *= e[(size_t)letter_mask(seq[i])].template segment<WINDOW>(len - i, WINDOW);

		if (pf >= p_mask) {
			if (mask_table)
				seq[i] = mask_table[(size_t)letter_mask(seq[i])];
			ranges.push_front(i);
		}

		t = f;
		t *= d;
		f *= f2f;
		f += p_repeat_end * b;
		b = b2b * b + t.sum();
	}
	return ranges;
}

}}}