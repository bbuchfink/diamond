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
A new repeat-masking method enables specific detection of homologous sequences, MC Frith, Nucleic Acids Research 2011 39(4):e23.*/

#include <array>
#include <stdint.h>
#include <algorithm>
#include <Eigen/Core>
#include "../basic/value.h"

using Eigen::Array;
using Eigen::Dynamic;

namespace Util { namespace tantan {

void mask(char *seq,
	int len,
	const float_t **likelihood_ratio_matrix,
	float_t p_repeat,
	float_t p_repeat_end,
	float_t p_mask,
	const char *mask_table) {
	constexpr int WINDOW = 50, RESERVE = 50000;

	thread_local std::array<Array<float_t, Dynamic, 1>, AMINO_ACID_COUNT> e;
	thread_local Array<float_t, Dynamic, 1> pb;
	thread_local Array<float_t, Dynamic, 1> scale;
	Array<float_t, WINDOW, 1> f(0.0);
	const float_t b2b = 1 - p_repeat, f2f = 1 - p_repeat_end;
	float_t b = 1.0;

	pb.resize(std::max(len, RESERVE));
	scale.resize(std::max((len - 1) / 16 + 1, (RESERVE - 1) / 16 + 1));

	for (int i = 0; i < AMINO_ACID_COUNT; ++i) {
		e[i].resize(std::max(RESERVE, len + WINDOW));
		const float_t *l = likelihood_ratio_matrix[i];
		float_t* p = &e[i][len - 1];
		for (int j = 0; j < len; ++j)
			*(p--) = l[seq[j]];
		std::fill(e[i].data() + len, e[i].data() + len + WINDOW, (float_t)0.0);
	}	

	for (int i = 0; i < len; ++i) {
		const float_t s = f.sum();
		f *= f2f;
		f += b * p_repeat / WINDOW;
		f *= e[(long)seq[i]].template segment<WINDOW>(len - i, WINDOW);
		b = b * b2b + s * p_repeat_end;

		if ((i & 15) == 15) {
			const float_t s = 1 / b;
			scale[i / 16] = s;
			b *= s;
			f *= s;
		}

		pb[i] = b;
	}

	const float_t z = b * b2b + f.sum() * p_repeat_end;
	b = b2b;
	f = p_repeat_end;

	for (int i = len - 1; i >= 0; --i) {
		const float_t pf = 1 - (pb[i] * b / z);

		if ((i & 15) == 15) {
			const float_t s = scale[i / 16];
			b *= s;
			f *= s;
		}

		f *= e[(long)seq[i]].template segment<WINDOW>(len - i, WINDOW);

		if (pf >= p_mask)
			seq[i] = mask_table[seq[i]];

		const float_t s = f.sum();
		f *= f2f;
		f += p_repeat_end * b;
		b = b2b * b + s * p_repeat / WINDOW;
	}
}

}}