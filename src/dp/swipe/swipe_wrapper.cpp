/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
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
#include "../../util/dynamic_iterator.h"

using std::list;
using std::atomic;
using std::thread;
using std::array;

namespace DP { namespace Swipe { namespace DISPATCH_ARCH {

template<typename _sv, typename _traceback, typename _cbs>
list<Hsp> swipe(const Sequence& query, Frame frame, DynamicIterator<DpTarget>& targets, _cbs composition_bias, vector<DpTarget>& overflow, Statistics& stats);

}}}

namespace DP { namespace BandedSwipe { namespace DISPATCH_ARCH {

template<typename _sv, typename _traceback, typename _cbs>
list<Hsp> swipe(
	const Sequence&query,
	Frame frame,
	vector<DpTarget>::const_iterator subject_begin,
	vector<DpTarget>::const_iterator subject_end,
	_cbs composition_bias,
	vector<DpTarget> &overflow,
	Statistics &stat);

template<typename _sv, typename _traceback>
list<Hsp> swipe_dispatch_cbs(
	const Sequence&query,
	Frame frame,
	vector<DpTarget>::const_iterator subject_begin,
	vector<DpTarget>::const_iterator subject_end,
	const int8_t* composition_bias,
	vector<DpTarget> &overflow,
	Statistics &stat)
{
	if (composition_bias == nullptr)
		return swipe<_sv, _traceback>(query, frame, subject_begin, subject_end, NoCBS(), overflow, stat);
	else
		return swipe<_sv, _traceback>(query, frame, subject_begin, subject_end, composition_bias, overflow, stat);
}

template<typename _sv, typename _traceback>
list<Hsp> full_swipe_dispatch_cbs(
	const Sequence&query,
	Frame frame,
	DynamicIterator<DpTarget>& targets,
	const int8_t* composition_bias,
	vector<DpTarget> &overflow,
	Statistics &stat)
{
	if (composition_bias == nullptr)
		return DP::Swipe::DISPATCH_ARCH::swipe<_sv, _traceback>(query, frame, targets, NoCBS(), overflow, stat);
	else
		return DP::Swipe::DISPATCH_ARCH::swipe<_sv, _traceback>(query, frame, targets, composition_bias, overflow, stat);
}

template<typename _sv>
list<Hsp> swipe_targets(const Sequence&query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	DynamicIterator<DpTarget>* targets,
	Frame frame,
	const int8_t *composition_bias,
	int flags,
	vector<DpTarget> &overflow,
	Statistics &stat)
{
	constexpr auto CHANNELS = vector<DpTarget>::const_iterator::difference_type(::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS);
	list<Hsp> out;
	if (flags & DP::FULL_MATRIX) {
		if (flags & TRACEBACK)
			return full_swipe_dispatch_cbs<_sv, VectorTraceback>(query, frame, *targets, composition_bias, overflow, stat);
		else if (flags & WITH_COORDINATES)
			return full_swipe_dispatch_cbs<_sv, ScoreWithCoords>(query, frame, *targets, composition_bias, overflow, stat);
		else
			return full_swipe_dispatch_cbs<_sv, ScoreOnly>(query, frame, *targets, composition_bias, overflow, stat);
	}
	else {
		for (vector<DpTarget>::const_iterator i = begin; i < end; i += CHANNELS) {
			if (flags & TRACEBACK) {
				if (config.traceback_mode == TracebackMode::VECTOR)
					out.splice(out.end(), swipe_dispatch_cbs<_sv, VectorTraceback>(query, frame, i, i + std::min(CHANNELS, end - i), composition_bias, overflow, stat));
				else
					out.splice(out.end(), swipe_dispatch_cbs<_sv, Traceback>(query, frame, i, i + std::min(CHANNELS, end - i), composition_bias, overflow, stat));
			}
			else
				out.splice(out.end(), swipe_dispatch_cbs<_sv, ScoreOnly>(query, frame, i, i + std::min(CHANNELS, end - i), composition_bias, overflow, stat));
		}
	}
	return out;
}

template<typename _sv>
void swipe_worker(const Sequence*query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	DynamicIterator<DpTarget>* targets,
	atomic<size_t> *next,
	Frame frame,
	const int8_t *composition_bias,
	int flags,
	list<Hsp> *out,
	vector<DpTarget> *overflow,
	Statistics *stat)
{
	Statistics stat2;
	size_t pos;
	vector<DpTarget> of;
	if (targets == nullptr) {
		while (begin + (pos = next->fetch_add(::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS)) < end)
			out->splice(out->end(), swipe_targets<_sv>(*query, begin + pos, std::min(begin + pos + ::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS, end), nullptr, frame, composition_bias, flags, of, stat2));
	}
	else
		out->splice(out->end(), swipe_targets<_sv>(*query, begin, end, targets, frame, composition_bias, flags, of, stat2));
	*overflow = std::move(of);
	*stat += stat2;
}

template<typename _sv>
list<Hsp> swipe_threads(const Sequence& query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	DynamicIterator<DpTarget>* targets,
	Frame frame,
	const int8_t *composition_bias,
	int flags,
	vector<DpTarget> &overflow,
	Statistics &stat) {
	if (end <= begin && (!targets || targets->count == 0))
		return {};

	std::unique_ptr<DynamicIterator<DpTarget>> my_targets;
	if (targets == nullptr && (flags & DP::FULL_MATRIX))
		my_targets.reset(new VectorIterator<DpTarget>(begin, end));

	if (flags & PARALLEL) {
		task_timer timer("Banded swipe (run)", config.target_parallel_verbosity);
		const size_t n = config.threads_align ? config.threads_align : config.threads_;
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
				targets ? targets : my_targets.get(),
				&next,
				frame,
				composition_bias,
				flags,
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
		return swipe_targets<_sv>(query, begin, end, targets ? targets : my_targets.get(), frame, composition_bias, flags, overflow, stat);
}

list<Hsp> recompute_reversed(const Sequence& query, Frame frame, const Bias_correction* composition_bias, int flags, Statistics& stat, list<Hsp>::const_iterator begin, list<Hsp>::const_iterator end) {
	array<vector<DpTarget>, 3> dp_targets;
	vector<DpTarget> overflow;
	const int qlen = (int)query.length();
#ifdef __SSE4_1__
	const int min_b = qlen <= UCHAR_MAX ? 0 : 1;
#else
	const int min_b = 1;
#endif
	for (auto i = begin; i != end; ++i) {
		int b = i->score <= UCHAR_MAX ? 0 : (i->score <= USHRT_MAX ? 1 : 2);
		b = std::max(b, min_b);
		dp_targets[b].emplace_back(i->target_seq, i->swipe_target);
	}

	list<Hsp> out;
#ifdef __SSE4_1__
	out = swipe_threads<::DISPATCH_ARCH::score_vector<int8_t>>(query, dp_targets[0].begin(), dp_targets[0].end(), nullptr, frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, overflow, stat);
#endif	
	out.splice(out.end(), swipe_threads<::DISPATCH_ARCH::score_vector<int16_t>>(query, dp_targets[1].begin(), dp_targets[1].end(), nullptr, frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, overflow, stat));
	out.splice(out.end(), swipe_threads<int32_t>(query, dp_targets[2].begin(), dp_targets[2].end(), nullptr, frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, overflow, stat));
	return out;
}

list<Hsp> swipe(const Sequence &query, vector<DpTarget> &targets8, const vector<DpTarget> &targets16, const vector<DpTarget>& targets32, DynamicIterator<DpTarget>* targets, Frame frame, const Bias_correction *composition_bias, int flags, Statistics &stat)
{
	vector<DpTarget> overflow8, overflow16, overflow32;
	list<Hsp> out;
	auto time_stat = (flags & TRACEBACK) ? Statistics::TIME_TRACEBACK_SW : Statistics::TIME_SW;
#ifdef __SSE4_1__
	if ((!targets8.empty() || targets) && config.cbs_matrix_scale < 16) {
		task_timer timer;
		std::sort(targets8.begin(), targets8.end());
		stat.inc(Statistics::TIME_TARGET_SORT, timer.microseconds());
		stat.inc(Statistics::EXT8, targets8.size());
		timer.go();
		out = swipe_threads<::DISPATCH_ARCH::score_vector<int8_t>>(query, targets8.begin(), targets8.end(), targets, frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, overflow8, stat);
		if ((flags & PARALLEL) == 0) stat.inc(time_stat, timer.microseconds());
	}
	else
		overflow8 = std::move(targets8);
#else
	overflow8 = std::move(targets8);
#endif
#ifdef __SSE2__
	if (!overflow8.empty() || !targets16.empty() || !targets32.empty()) {
		overflow8.insert(overflow8.end(), targets16.begin(), targets16.end());
		stat.inc(Statistics::EXT16, overflow8.size());
		task_timer timer;
		std::sort(overflow8.begin(), overflow8.end());
		stat.inc(Statistics::TIME_TARGET_SORT, timer.microseconds());
		timer.go();
		out.splice(out.end(), swipe_threads<::DISPATCH_ARCH::score_vector<int16_t>>(query, overflow8.begin(), overflow8.end(), nullptr, frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, overflow16, stat));
		if ((flags & PARALLEL) == 0) stat.inc(time_stat, timer.microseconds());
		if (!overflow16.empty() || !targets32.empty()) {
			overflow16.insert(overflow16.end(), targets32.begin(), targets32.end());
			stat.inc(Statistics::EXT32, overflow16.size());
			timer.go();
			out.splice(out.end(), swipe_threads<int32_t>(query, overflow16.begin(), overflow16.end(), nullptr, frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, overflow32, stat));
			stat.inc(time_stat, timer.microseconds());
		}
	}
	return out;
#else
	overflow8.insert(overflow8.end(), targets16.begin(), targets16.end());
	stat.inc(Statistics::EXT32, overflow8.size());
	return swipe_threads<int32_t>(query, overflow8.begin(), overflow8.end(), targets, frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, overflow32, stat);
#endif
}
		
}}}
