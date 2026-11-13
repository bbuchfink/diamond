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
#include <chrono>
#include <mutex>
#include <string>
#include <vector>
#include <unordered_set>
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

/* A list of integers stored in a file, one number per line. The file is only ever appended to,
   so the numbers already read are cached and each update only reads the new part of the file. */
struct AtomicList {

	template<typename Job>
	AtomicList(const std::string& file_name, Job& job):
		stack_(file_name, job),
		offset_(0)
	{
	}

	void push(int64_t i) {
		stack_.push(std::to_string(i) + '\n');
	}

	/* Waits until the list holds at least n numbers. Returns false without waiting further if the job
	   has been stopped (Job::stopped()); the caller then has to quit. */
	template<typename Job>
	bool await(int64_t n, Job& job) {
		for (;;) {
			if (update() >= n)
				return true;
			if (job.stopped())
				return false;
			std::this_thread::sleep_for(std::chrono::seconds(1));
		}
	}

	// Returns the numbers currently in the list.
	std::vector<int64_t> values() {
		update();
		std::lock_guard<std::mutex> lock(mtx_);
		return values_;
	}

	void close() {
		stack_.close();
	}

	void remove() {
		stack_.remove();
	}

private:

	// Reads the lines appended since the last update and returns the number of values.
	int64_t update() {
		std::lock_guard<std::mutex> lock(mtx_);
		const std::string buf = stack_.read_from(offset_);
		size_t begin = 0, end;
		while ((end = buf.find('\n', begin)) != std::string::npos) {
			values_.push_back(std::stoll(buf.substr(begin, end - begin)));
			begin = end + 1;
		}
		offset_ += begin;
		return (int64_t)values_.size();
	}

	FileStack stack_;
	int64_t offset_;
	std::vector<int64_t> values_;
	std::mutex mtx_;

};

struct Atomic {

	template<typename Job>
	Atomic(const std::string& file_name, Job& job):
		stack_(file_name, job)
	{
	}

	/* The numbers in skip are read once here; fetch_add skips over them. Numbers added to the
	   list later are not tracked. */
	template<typename Job>
	Atomic(const std::string& file_name, Job& job, AtomicList& skip):
		stack_(file_name, job)
	{
		const std::vector<int64_t> v = skip.values();
		skip_.insert(v.begin(), v.end());
	}

	int64_t get() {
		int64_t i;
		stack_.top(i);
		return i >= 0 ? i : 0;
	}

	int64_t fetch_add(int64_t n = 1) {
		int64_t i;
		do
			i = stack_.fetch_add(n);
		while (skip_.find(i) != skip_.end());
		return i;
	}

	/* Waits until the value is at least n. Returns false without waiting further if the job
	   has been stopped (Job::stopped()); the caller then has to quit. */
	template<typename Job>
	bool await(int64_t n, Job& job) {
		for (;;) {
			int64_t i;
			stack_.top(i);
			if (i >= n)
				return true;
			if (job.stopped())
				return false;
			std::this_thread::sleep_for(std::chrono::seconds(1));
		}
	}

	// Resets the value to 0.
	void reset() {
		stack_.clear();
	}

	void close() {
		stack_.close();
	}

	void remove() {
		stack_.remove();
	}

private:

	FileStack stack_;
	std::unordered_set<int64_t> skip_;

};