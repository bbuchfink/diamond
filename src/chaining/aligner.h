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
#include <map>
#include "diag_graph.h"
#include "basic/match.h"

namespace Chaining {

struct Aligner
{

	enum { link_padding = 10, reverse_link_min_overhang = 10 };

	int get_approximate_link(int d_idx, int e_idx, double space_penalty, int max_i);
	template<typename _it>
	void forward_pass(_it begin, _it end, bool init, double space_penalty);
	void backtrace(const size_t node, const int j_end, Hsp* out, ApproxHsp& t, const int score_max, const int score_min, const int max_shift, unsigned& next) const;
	bool backtrace_old(size_t node, int j_end, Hsp* out, ApproxHsp& t, int score_max, int score_min, int max_shift, unsigned& next) const;
	void backtrace(size_t top_node, Hsp* out, ApproxHsp& t, int max_shift, unsigned& next, int max_j) const;
	int backtrace(size_t top_node, std::list<Hsp>& hsps, std::list<ApproxHsp>& ts, std::list<ApproxHsp>::iterator& t_begin, int cutoff, int max_shift) const;
	int backtrace(std::list<Hsp>& hsps, std::list<ApproxHsp>& ts, int cutoff, int max_shift) const;
	int run(std::list<Hsp>& hsps, std::list<ApproxHsp>& ts, double space_penalty, int cutoff, int max_shift);
	int run(std::list<Hsp>& hsps, std::list<ApproxHsp>& ts, std::vector<DiagonalSegment>::const_iterator begin, std::vector<DiagonalSegment>::const_iterator end, int band);
	Aligner(const Sequence& query, const Sequence& subject, bool log, unsigned frame);

	const Sequence query, subject;
	//const Bias_correction &query_bc;
	const bool log;
	const unsigned frame;
	static thread_local DiagGraph diags;
	static thread_local std::map<int, unsigned> window;

};

}