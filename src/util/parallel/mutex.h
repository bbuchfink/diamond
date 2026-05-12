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
#include <mutex>

struct Sync {};
struct Async {};

template<typename Tag>
struct Mutex {};

template<>
struct Mutex<Sync> {
	void lock() {}
	void unlock() {}
};

template<>
struct Mutex<Async> {
	void lock() {
		mtx_.lock();
	}
	void unlock() {
		mtx_.unlock();
	}
private:
	std::mutex mtx_;
};