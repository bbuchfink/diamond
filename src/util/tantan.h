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

#ifndef TANTAN_H_
#define TANTAN_H_

#include "simd.h"

namespace Util { namespace tantan {

typedef float float_t;

DECL_DISPATCH(void, mask, (char *seq, int len, const float_t **likelihood_ratio_matrix, float_t p_repeat, float_t p_repeat_end, float_t repeat_decay, float_t p_mask, const char *maskTable))

}}

#endif