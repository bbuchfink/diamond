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
#include <mutex>
#include <thread>
#if defined(__SSE2__) || defined(_M_AMD64)
#include <emmintrin.h>
#endif

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

// Lock for short critical sections, which take less time than suspending and waking up a thread.
// A waiting thread spins instead of sleeping. It yields to the scheduler if the lock is not released soon, as the
// holder may have been preempted.
struct SpinLock {
	void lock() noexcept {
		int spins = 0;
		while (locked_.exchange(true, std::memory_order_acquire))
			while (locked_.load(std::memory_order_relaxed))
				if (++spins < SPINS_BEFORE_YIELD)
					pause();
				else
					std::this_thread::yield();
	}
	bool try_lock() noexcept {
		return !locked_.load(std::memory_order_relaxed) && !locked_.exchange(true, std::memory_order_acquire);
	}
	void unlock() noexcept {
		locked_.store(false, std::memory_order_release);
	}
private:
	static constexpr int SPINS_BEFORE_YIELD = 256;
	static void pause() noexcept {
#if defined(__SSE2__) || defined(_M_AMD64)
		_mm_pause();
#elif defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
		__asm__ __volatile__("yield");
#endif
	}
	std::atomic<bool> locked_{ false };
};
