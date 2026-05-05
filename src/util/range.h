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

template<typename It>
struct Range {

	Range() {}

	Range(const It& begin, const It& end):
		begin_(begin),
		end_(end)
	{}

	It begin() const {
		return begin_;
	}

	It end() const {
		return end_;
	}

	size_t size() const {
		return end_ - begin_;
	}

private:

	It begin_, end_;

};