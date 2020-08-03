/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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


#include <list>
#include <atomic>
#include <thread>
#include <numeric>
#include <limits.h>
#include "../dp.h"
#include "../score_vector_int16.h"
#include "../score_vector_int8.h"
#include "../../util/log_stream.h"

using std::list;
using std::atomic;
using std::thread;

namespace DP { namespace Swipe { namespace DISPATCH_ARCH {

template<typename _sv, typename _traceback, typename _cbs>
list<Hsp> swipe(const sequence& query, Frame frame, vector<DpTarget>::const_iterator subject_begin, vector<DpTarget>::const_iterator subject_end, _cbs composition_bias, int score_cutoff, vector<DpTarget>& overflow, Statistics& stats);

}}}

namespace DP { namespace BandedSwipe { namespace DISPATCH_ARCH {

template<typename _sv, typename _traceback, typename _cbs>
list<Hsp> swipe(
	const sequence &query,
	Frame frame,
	vector<DpTarget>::const_iterator subject_begin,
	vector<DpTarget>::const_iterator subject_end,
	_cbs composition_bias,
	int score_cutoff,
	vector<DpTarget> &overflow,
	Statistics &stat);

template<typename _sv>
list<Hsp> swipe_targets(const sequence &query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	Frame frame,
	const int8_t *composition_bias,
	int flags,
	int score_cutoff,
	vector<DpTarget> &overflow,
	Statistics &stat)
{
	constexpr auto CHANNELS = vector<DpTarget>::const_iterator::difference_type(::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS);
	list<Hsp> out;
	if (flags & DP::FULL_MATRIX) {
		return DP::Swipe::DISPATCH_ARCH::swipe<_sv, ScoreOnly>(query, frame, begin, end, NoCBS(), score_cutoff, overflow, stat);
	}
	else {
		for (vector<DpTarget>::const_iterator i = begin; i < end; i += CHANNELS) {
			if (flags & TRACEBACK) {
				if (config.traceback_mode == TracebackMode::STAT) {
					if (composition_bias == nullptr)
						out.splice(out.end(), swipe<_sv, StatTraceback>(query, frame, i, i + std::min(CHANNELS, end - i), NoCBS(), score_cutoff, overflow, stat));
					else
						out.splice(out.end(), swipe<_sv, StatTraceback>(query, frame, i, i + std::min(CHANNELS, end - i), composition_bias, score_cutoff, overflow, stat));
				}
				else if (config.traceback_mode == TracebackMode::VECTOR) {
					if (composition_bias == nullptr)
						out.splice(out.end(), swipe<_sv, VectorTraceback>(query, frame, i, i + std::min(CHANNELS, end - i), NoCBS(), score_cutoff, overflow, stat));
					else
						out.splice(out.end(), swipe<_sv, VectorTraceback>(query, frame, i, i + std::min(CHANNELS, end - i), composition_bias, score_cutoff, overflow, stat));
				}
				else {
					if (composition_bias == nullptr)
						out.splice(out.end(), swipe<_sv, Traceback>(query, frame, i, i + std::min(CHANNELS, end - i), NoCBS(), score_cutoff, overflow, stat));
					else
						out.splice(out.end(), swipe<_sv, Traceback>(query, frame, i, i + std::min(CHANNELS, end - i), composition_bias, score_cutoff, overflow, stat));
				}
			}
			else {
				if (composition_bias == nullptr)
					out.splice(out.end(), swipe<_sv, ScoreOnly>(query, frame, i, i + std::min(CHANNELS, end - i), NoCBS(), score_cutoff, overflow, stat));
				else
					out.splice(out.end(), swipe<_sv, ScoreOnly>(query, frame, i, i + std::min(CHANNELS, end - i), composition_bias, score_cutoff, overflow, stat));
			}
		}
	}
	return out;
}

template<typename _sv>
void swipe_worker(const sequence *query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	atomic<size_t> *next,
	Frame frame,
	const int8_t *composition_bias,
	int flags,
	int score_cutoff,
	list<Hsp> *out,
	vector<DpTarget> *overflow,
	Statistics *stat)
{
	Statistics stat2;
	size_t pos;
	vector<DpTarget> of;
	while (begin + (pos = next->fetch_add(::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS)) < end)
		out->splice(out->end(), swipe_targets<_sv>(*query, begin + pos, std::min(begin + pos + ::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS, end), frame, composition_bias, flags, score_cutoff, of, stat2));
	*overflow = std::move(of);
	*stat += stat2;
}

template<typename _sv>
list<Hsp> swipe_threads(const sequence &query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	Frame frame,
	const int8_t *composition_bias,
	int flags,
	int score_cutoff,
	vector<DpTarget> &overflow,
	Statistics &stat) {
	if (flags & PARALLEL) {
		task_timer timer("Banded swipe (run)", config.target_parallel_verbosity);
		const size_t n = config.threads_;
		vector<thread> threads;
		vector<list<Hsp>> thread_out(n);
		vector<vector<DpTarget>> thread_overflow(n);
		atomic<size_t> next(0);
		for (size_t i = 0; i < n; ++i)
			threads.emplace_back(
				swipe_worker<_sv>,
				&query,
				begin,
				end,
				&next,
				frame,
				composition_bias,
				flags,
				score_cutoff,
				&thread_out[i],
				&thread_overflow[i],
				&stat);
		for (auto &t : threads)
			t.join();
		timer.go("Banded swipe (merge)");
		list<Hsp> out;
		for (list<Hsp> &l : thread_out)
			out.splice(out.end(), l);
		overflow.reserve(std::accumulate(thread_overflow.begin(), thread_overflow.end(), (size_t)0, [](size_t n, const vector<DpTarget> &v) { return n + v.size(); }));
		for (const vector<DpTarget> &v : thread_overflow)
			overflow.insert(overflow.end(), v.begin(), v.end());
		return out;
	}
	else
		return swipe_targets<_sv>(query, begin, end, frame, composition_bias, flags, score_cutoff, overflow, stat);
}

list<Hsp> swipe(const sequence &query, vector<DpTarget> &targets8, vector<DpTarget> &targets16, Frame frame, const Bias_correction *composition_bias, int flags, int score_cutoff, Statistics &stat)
{
	vector<DpTarget> overflow8, overflow16, overflow32;
	list<Hsp> out;
	auto time_stat = (flags & TRACEBACK) ? Statistics::TIME_TRACEBACK_SW : Statistics::TIME_SW;
#ifdef __SSE4_1__
	task_timer timer;
	std::sort(targets8.begin(), targets8.end());
	stat.inc(Statistics::TIME_TARGET_SORT, timer.microseconds());
	stat.inc(Statistics::EXT8, targets8.size());
	timer.go();
	out = swipe_threads<::DISPATCH_ARCH::score_vector<int8_t>>(query, targets8.begin(), targets8.end(), frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, score_cutoff, overflow8, stat);
	if((flags & PARALLEL) == 0) stat.inc(time_stat, timer.microseconds());
#else
	overflow8 = std::move(targets8);
#endif
#ifdef __SSE2__
	if (!overflow8.empty() || !targets16.empty()) {
		overflow8.insert(overflow8.end(), targets16.begin(), targets16.end());
		stat.inc(Statistics::EXT16, overflow8.size());
		task_timer timer;
		std::sort(overflow8.begin(), overflow8.end());
		stat.inc(Statistics::TIME_TARGET_SORT, timer.microseconds());
		timer.go();
		out.splice(out.end(), swipe_threads<::DISPATCH_ARCH::score_vector<int16_t>>(query, overflow8.begin(), overflow8.end(), frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, score_cutoff, overflow16, stat));
		if ((flags & PARALLEL) == 0) stat.inc(time_stat, timer.microseconds());
		if (!overflow16.empty()) {
			stat.inc(Statistics::EXT32, overflow16.size());
			timer.go();
			out.splice(out.end(), swipe_threads<int32_t>(query, overflow16.begin(), overflow16.end(), frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, score_cutoff, overflow32, stat));
			stat.inc(time_stat, timer.microseconds());
		}
	}
	return out;
#else
	overflow8.insert(overflow8.end(), targets16.begin(), targets16.end());
	stat.inc(Statistics::EXT32, overflow8.size());
	return swipe_threads<int32_t>(query, overflow8.begin(), overflow8.end(), frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, score_cutoff, overflow32, stat);
#endif
}
		
}}}
