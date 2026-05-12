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
#include <atomic>
#include <thread>
#include "filestack.h"

template<typename T>
void exchange_if_smaller(std::atomic<T>& value, T replacement) {
	T current = value.load(std::memory_order_relaxed);
	while (current > replacement &&
		!value.compare_exchange_weak(
			current,
			replacement,
			std::memory_order_relaxed,
			std::memory_order_relaxed)) {
	}
}

template<typename T>
void exchange_if_larger(std::atomic<T>& value, T replacement) {
	T current = value.load(std::memory_order_relaxed);
	while (current < replacement &&
		!value.compare_exchange_weak(
			current,
			replacement,
			std::memory_order_relaxed,
			std::memory_order_relaxed)) {
	}
}

struct Atomic {

	template<typename Job>
	Atomic(const std::string& file_name, Job& job):
		stack_(file_name, job)
	{
	}

	int64_t get() {
		int64_t i;
		stack_.top(i);
		return i >= 0 ? i : 0;
	}

	int64_t fetch_add(int64_t n = 1) {
		return stack_.fetch_add(n);
	}

	void await(int64_t n) {
		for (;;) {
			int64_t i;
			stack_.top(i);
			if (i >= n)
				return;
			std::this_thread::sleep_for(std::chrono::seconds(1));
		}
	}

	void close() {
		stack_.close();
	}

	void remove() {
		stack_.remove();
	}

private:

	FileStack stack_;

};