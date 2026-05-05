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
#include "log_stream.h"

struct Profiler {

	Profiler(const char* key) :
		timer(),
		key(key)
	{}

	~Profiler() {
		finish();
	}

	void finish() {
		if (!key) return;
		times[key] += timer.nanoseconds();
		key = nullptr;
	}

	static void print(size_t n) {
		for (const auto& i : times)
			message_stream << i.first << ": " << (double)i.second/n/1e3 << " micros" << std::endl;
	}

	TaskTimer timer;
	const char* key;
	static std::map<std::string, uint64_t> times;

};