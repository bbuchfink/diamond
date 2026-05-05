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

#pragma once
#include "basic/match.h"
#include "output/output_format.h"

namespace Extension {

bool filter_hsp(Hsp& hsp,
	int source_query_len,
	const char *query_title,
	int subject_len,
	const char* subject_title,
	const Sequence& query_seq,
	const Sequence& subject_seq,
	const double query_self_aln_score,
	const double target_self_aln_score,
	const OutputFormat* output_format);

}