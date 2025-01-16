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
#include <mutex>
#include <numeric>
#include "../score_vector_int16.h"
#include "../score_vector_int8.h"
#include "util/simd/vector.h"
#include "../dp.h"
#include "util/log_stream.h"
#include "data/sequence_set.h"
#include "../score_vector.h"
#include "cell_update.h"
#include "full_swipe.h"
#include "util/geo/geo.h"
#include "util/simd/dispatch.h"
#include "banded_swipe.h"
#include "stats/hauser_correction.h"
#include "stats/stats.h"

using std::list;
using std::atomic;
using std::thread;
using std::array;
using std::atomic_size_t;
using std::mutex;
using std::accumulate;
using std::string;
using std::pair;

template<bool tb, typename RC, typename C, typename IdM>
struct SwipeConfig {
	static constexpr bool traceback = tb;
	using RowCounter = RC;
	using Cell = C;
	using IdMask = IdM;
};

namespace DP { namespace BandedSwipe { namespace DISPATCH_ARCH {

static void sort(const TargetVec::iterator begin, const TargetVec::iterator end) {
	std::sort(begin, end);
}

static void sort(const SequenceSet::ConstIterator begin, const SequenceSet::ConstIterator end) {
}

static unsigned bin(int x) {
	return x < UCHAR_MAX ? 0 : (x < USHRT_MAX ? 1 : 2);
}

static const HspValues NO_TRACEBACK = HspValues::COORDS | HspValues::IDENT | HspValues::LENGTH | HspValues::MISMATCHES | HspValues::GAP_OPENINGS;

unsigned bin(HspValues v, int query_len, int score, int ungapped_score, const int64_t dp_size, unsigned score_width, const Loc mismatch_est) {
	unsigned b = 0;
	b = std::max(b, bin(score));
	if (ungapped_score > config.cutoff_score_8bit)
		b = std::max(b, 1u);
	b = std::max(b, score_width);
	b = std::max(b, bin(mismatch_est));
#ifndef __ARM_NEON
#ifndef __SSE4_1__
	b = std::max(b, 1u);
#endif
#ifndef __SSE2__
	b = 2;
#endif
#endif
	if (v != HspValues::NONE) {
		b = std::max(b, bin(query_len));
		if (dp_size > config.max_swipe_dp) {
			if (flag_only(v, NO_TRACEBACK))
				b += SCORE_BINS;
			else
				b = 2;
		}
		else if (flag_only(v, HspValues::COORDS) && !config.approx_backtrace)
			b += SCORE_BINS;
	}
	return b;
}

template<typename Sv>
static int64_t matrix_size(const int query_len, const TargetVec::const_iterator begin, const TargetVec::const_iterator end, const Flags flags) {
	int64_t s = 0;
	for (auto i = begin; i != end; ++i) {
		const int64_t cols = flag_any(flags, Flags::FULL_MATRIX) ? i->seq.length() : i->cols,
			size = int64_t(flag_any(flags, Flags::FULL_MATRIX) ? query_len : i->d_end - i->d_begin) * cols * ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS / 2;
		s = std::max(size, s);
	}
	return s;
}

template<typename Sv>
static size_t matrix_size(const int query_len, const SequenceSet::ConstIterator begin, const SequenceSet::ConstIterator end, const Flags flags) {
	return 0;
}

static bool reversed(const HspValues v) {
	return flag_only(v, NO_TRACEBACK)
		&& flag_any(v, HspValues::QUERY_START | HspValues::TARGET_START | HspValues::MISMATCHES | HspValues::GAP_OPENINGS);
}

template<typename Sv, typename Cbs, typename Cfg>
static list<Hsp> dispatch_swipe(const TargetVec::const_iterator subject_begin, const TargetVec::const_iterator subject_end, Cbs composition_bias, TargetVec& overflow, Params& p)
{
	return ::DP::BandedSwipe::DISPATCH_ARCH::swipe<Sv, Cbs, Cfg>(subject_begin, subject_end, composition_bias, overflow, p);
}

template<typename Sv, typename Cbs, typename Cfg>
static list<Hsp> dispatch_swipe(const SequenceSet::ConstIterator subject_begin, const SequenceSet::ConstIterator subject_end, Cbs composition_bias, TargetVec& overflow, Params& p)
{
	return {};
}

template<typename Sv, typename Cbs, typename It, typename Cfg>
static list<Hsp> dispatch_swipe(const It begin, const It end, atomic<BlockId>* const next, Cbs composition_bias, TargetVec& overflow, Params& p)
{
	constexpr auto CHANNELS = TargetVec::const_iterator::difference_type(::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS);
	if (flag_any(p.flags, Flags::FULL_MATRIX))
		return ::DP::Swipe::DISPATCH_ARCH::swipe<Sv, Cbs, It, Cfg>(begin, end, next, composition_bias, overflow, p);
	else {
		list<Hsp> out;
		for (It i = begin; i < end; i += std::min(CHANNELS, end - i))
			out.splice(out.end(), dispatch_swipe<Sv, Cbs, Cfg>(i, i + std::min(CHANNELS, end - i), composition_bias, overflow, p));
		return out;
	}
}

template<typename Sv, typename It, typename Cfg>
static list<Hsp> dispatch_swipe(const It begin, const It end, atomic<BlockId>* const next, TargetVec &overflow, Params& p)
{
	if (p.composition_bias == nullptr)
		return dispatch_swipe<Sv, NoCBS, It, Cfg>(begin, end, next, NoCBS(), overflow, p);
	else
		return dispatch_swipe<Sv, const int8_t*, It, Cfg>(begin, end, next, p.composition_bias, overflow, p);
}

template<typename Sv, typename It>
static list<Hsp> dispatch_swipe(const It begin, const It end, atomic<BlockId>* const next, TargetVec &overflow, const int round, const int bin, Params& p)
{
	if (p.v == HspValues::NONE) {
		using Cfg = SwipeConfig<false, DummyRowCounter<Sv>, Sv, DummyIdMask<Sv>>;
		return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
	}
	if (bin < SCORE_BINS) {
		using Cfg = SwipeConfig<true, VectorRowCounter<Sv>, Sv, DummyIdMask<Sv>>;
		return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
	}
	if (round == 0) {
		if (!flag_any(p.v, HspValues::IDENT | HspValues::LENGTH)) {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, Sv, DummyIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
		}
		else {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, ForwardCell<Sv>, VectorIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
		}
	}
	else if (round == 1) {
		if (!flag_any(p.v, HspValues::MISMATCHES | HspValues::GAP_OPENINGS)) {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, Sv, DummyIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
		}
		else {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, BackwardCell<Sv>, VectorIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
		}
	}
	throw std::runtime_error("Unreachable");
}

template<typename Sv, typename It>
static void swipe_worker(const It begin, const It end, atomic<BlockId>* const next, list<Hsp> *out, TargetVec *overflow, const int round, const int bin, Params* p)
{
	const ptrdiff_t CHANNELS = ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS;
	Statistics stat2;
	size_t pos;
	TargetVec of;
	Params params{
		p->query,
		p->query_id,
		p->frame,
		p->query_source_len,
		p->composition_bias,
		p->flags,
		p->reverse_targets,
		p->target_max_len,
		p->v,
		stat2,
		nullptr
	};
	if (flag_any(p->flags, Flags::FULL_MATRIX))
		*out = dispatch_swipe<Sv, It>(begin, end, next, of, round, bin, params);
	else
		while (begin + (pos = next->fetch_add(CHANNELS, std::memory_order_relaxed)) < end) {
			const auto start = begin + pos;
			out->splice(out->end(), dispatch_swipe<Sv, It>(start, start + std::min(CHANNELS, end - start), next, of, round, bin, params));
		}
		
	*overflow = std::move(of);
	p->stat += stat2;
}

template<typename Sv, typename It>
static void swipe_task(const It begin, const It end, list<Hsp> *out, TargetVec *overflow, mutex* mtx, const int round, const int bin, Params* p) {
	const ptrdiff_t CHANNELS = ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS;
	Statistics stat2;
	TargetVec of;
	atomic<BlockId> next(0);
	Params params{
		p->query,
		p->query_id,
		p->frame,
		p->query_source_len,
		p->composition_bias,
		p->flags,
		p->reverse_targets,
		p->target_max_len,
		p->v,
		stat2,
		nullptr
	};
	list<Hsp> hsp = dispatch_swipe<Sv, It>(begin, end, &next, of, round, bin, params);
	{
		std::lock_guard<mutex> lock(*mtx);
		overflow->push_back(of);
		out->splice(out->end(), hsp);
	}
	p->stat += stat2;
}

template<typename Sv, typename It>
static list<Hsp> swipe_threads(const It begin, const It end, TargetVec &overflow, const int round, const int bin, Params& p) {
	const ptrdiff_t CHANNELS = ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS;
	if (begin == end)
		return {};

	atomic<BlockId> next(0);
	if (flag_any(p.flags, Flags::PARALLEL)) {
		TaskTimer timer("Banded swipe (run)", config.target_parallel_verbosity);
		const size_t n = config.threads_align ? config.threads_align : config.threads_;
		vector<thread> threads;
		vector<list<Hsp>> thread_out(n);
		vector<TargetVec> thread_overflow(n);
		for (size_t i = 0; i < n; ++i)
			threads.emplace_back(swipe_worker<Sv, It>, begin, end, &next, &thread_out[i], &thread_overflow[i], round, bin, &p);
		for (auto &t : threads)
			t.join();
		timer.go("Banded swipe (merge)");
		list<Hsp> out;
		for (list<Hsp> &l : thread_out)
			out.splice(out.end(), l);
		overflow.reserve(std::accumulate(thread_overflow.begin(), thread_overflow.end(), (size_t)0, [](size_t n, const TargetVec &v) { return n + v.size(); }));
		for (const TargetVec& v : thread_overflow)
			overflow.push_back(v);
		return out;
	}

	if(!p.thread_pool)
		return dispatch_swipe<Sv, It>(begin, end, &next, overflow, round, bin, p);

	list<Hsp> hsp;
	ThreadPool::TaskSet task_set(*p.thread_pool, 0);
	mutex mtx;
	int64_t size = 0;
	It i0 = begin, i1 = begin;
	while(i1 < end) {
		const auto n = std::min(CHANNELS, end - i1);
		size += accumulate(i1, i1 + n, (int64_t)0, [&p](int64_t n, const DpTarget& t) {return n + t.cells(p.flags, p.query.length()); });
		i1 += n;
		if (size >= config.swipe_task_size) {
			task_set.enqueue(swipe_task<Sv, It>, i0, i1, &hsp, &overflow, &mtx, round, bin, &p);
			p.stat.inc(Statistics::SWIPE_TASKS_TOTAL);
			p.stat.inc(Statistics::SWIPE_TASKS_ASYNC);
			i0 = i1;
			size = 0;
		}
	}
	if (task_set.total() == 0) {
		p.stat.inc(Statistics::SWIPE_TASKS_TOTAL);
		return dispatch_swipe<Sv, It>(i0, i1, &next, overflow, round, bin, p);
	}
	if (i1 - i0 > 0) {
		p.stat.inc(Statistics::SWIPE_TASKS_TOTAL);
		p.stat.inc(Statistics::SWIPE_TASKS_ASYNC);
		task_set.enqueue(swipe_task<Sv, It>, i0, i1, &hsp, &overflow, &mtx, round, bin, &p);
	}
	task_set.run();
	return hsp;
}

template<typename It>
static pair<list<Hsp>, TargetVec> swipe_bin(const unsigned bin, const It begin, const It end, const int round, Params& p) {
	if (end - begin == 0)
		return { {},{} };
	TargetVec overflow;
	list<Hsp> out;
	auto time_stat = flag_any(p.v, HspValues::TRANSCRIPT) ? Statistics::TIME_TRACEBACK_SW : Statistics::TIME_SW;
	if (!flag_any(p.flags, Flags::FULL_MATRIX))
		sort(begin, end);
	p.stat.inc(Statistics::value(Statistics::EXT8 + bin % SCORE_BINS), end - begin);
	TaskTimer timer;
	switch (bin) {
#if defined(__SSE4_1__) | defined(__ARM_NEON)
	case 0:
	case 3:
		if (flag_any(p.flags, Flags::SEMI_GLOBAL))
			out = swipe_threads<::DISPATCH_ARCH::ScoreVector<int8_t, 0>, It>(begin, end, overflow, round, bin, p);
		else
			out = swipe_threads<::DISPATCH_ARCH::ScoreVector<int8_t, SCHAR_MIN>, It>(begin, end, overflow, round, bin, p);
		break;
#endif
#if defined(__SSE2__) | defined(__ARM_NEON)
	case 1:
	case 4:
		if (flag_any(p.flags, Flags::SEMI_GLOBAL))
			out = swipe_threads<::DISPATCH_ARCH::ScoreVector<int16_t, 0>, It>(begin, end, overflow, round, bin, p);
		else
			out = swipe_threads<::DISPATCH_ARCH::ScoreVector<int16_t, SHRT_MIN>, It>(begin, end, overflow, round, bin, p);
		break;
#endif
	case 2:
	case 5:
		out = swipe_threads<int32_t, It>(begin, end, overflow, round, bin, p);
		break;
	default:
		throw std::runtime_error("Invalid SWIPE bin.");
	}
	if (!flag_any(p.flags, Flags::PARALLEL)) p.stat.inc(time_stat, timer.microseconds());
	return { out, overflow };
}

static Loc mismatch_est(const Loc query_len, const Loc target_len, const int32_t aln_len, const HspValues v) {
	if (!flag_any(v, HspValues::MISMATCHES))
		return 0;
	const Loc m = std::min(query_len, target_len);
	return aln_len > 0 ? std::min(aln_len, m) : m;
}

static list<Hsp> recompute_reversed(list<Hsp> &hsps, Params& p) {
	Targets dp_targets;
	const Loc qlen = p.query.length();

	for (auto i = hsps.begin(); i != hsps.end(); ++i) {
		const double qcov = i->query_cover_percent(p.query_source_len), tcov = i->subject_cover_percent(i->target_seq.length());
		const pair<Loc, Loc> min_range_len = i->min_range_len(config.query_or_target_cover > 0 ? config.query_or_target_cover : config.query_cover,
			config.query_or_target_cover > 0 ? config.query_or_target_cover : config.subject_cover, qlen, i->target_seq.length());
		const double qa = Stats::approx_id(i->score, min_range_len.first, 0), ta = Stats::approx_id(i->score, min_range_len.second, 0);
		if (qcov < config.query_cover || tcov < config.subject_cover || std::max(qcov, tcov) < config.query_or_target_cover
			|| (config.query_or_target_cover == 0 && std::min(qa, ta) < config.approx_min_id.get(0.0))
			|| (config.query_or_target_cover > 0 && std::max(qa, ta) < config.approx_min_id.get(0.0)))
			continue;
		Sequence reversed_seq(i->target_seq.data(), i->target_seq.data() + i->subject_range.end_);
		const Loc band = flag_any(p.flags, Flags::FULL_MATRIX) ? qlen : i->d_end - i->d_begin,
			tlen = i->subject_range.end_,
			b = bin(p.v, band, i->score, 0, INT64_MAX, 0, mismatch_est(i->query_range.end_, tlen, i->length, p.v));
		assert(b >= SCORE_BINS);
		const DpTarget::CarryOver carry_over{ i->query_range.end_, i->subject_range.end_, i->identities, i->length };
		dp_targets[b].emplace_back(reversed_seq, i->target_seq.length(), Geo::rev_diag(i->d_end-1, qlen, tlen), Geo::rev_diag(i->d_begin, qlen, tlen) + 1,
			Interval(), 0, i->swipe_target, qlen, i->matrix, carry_over);
	}

	vector<Letter> reversed = p.query.reverse();
	vector<int8_t> rev_cbs = HauserCorrection::reverse(p.composition_bias, p.query.length());
	const int8_t* cbs = p.composition_bias ? rev_cbs.data() : nullptr;
	Params params{
		Sequence(reversed),
		p.query_id,
		p.frame,
		p.query_source_len,
		cbs,
		p.flags,
		true,
		0,
		p.v,
		p.stat,
		p.thread_pool
	};
	list<Hsp> out;
	for (unsigned bin = SCORE_BINS; bin < BINS; ++bin) {
		params.target_max_len = dp_targets[bin].max_len();
		auto r = swipe_bin(bin, dp_targets[bin].begin(), dp_targets[bin].end(), 1, params);
		if (!r.second.empty())
			throw std::runtime_error("Non-empty overflow list in reversed DP. Query = " + string(p.query_id) + " bin=" + std::to_string(bin)
				+ " target=" + r.second.front().seq.to_string()
				+ " d_begin=" + std::to_string(r.second.front().d_begin)
				+ " d_end=" + std::to_string(r.second.front().d_end));
		out.splice(out.end(), r.first);
	}
	return out;
}

list<Hsp> swipe(const Targets &targets, Params& p)
{
	pair<list<Hsp>, TargetVec> result;
	list<Hsp> out, out_tmp;
	for (int algo_bin = 0; algo_bin < ALGO_BINS; ++algo_bin) {
		for (int score_bin = 0; score_bin < SCORE_BINS; ++score_bin) {
			const int bin = algo_bin * SCORE_BINS + score_bin;
			TargetVec round_targets;
			round_targets.reserve(targets[bin].size() + result.second.size());
			round_targets.push_back(targets[bin]);
			round_targets.push_back(result.second);
			p.target_max_len = round_targets.max_len();
			result = swipe_bin(bin, round_targets.begin(), round_targets.end(), 0, p);
			if(algo_bin == 0)
				out.splice(out.end(), result.first);
			else
				out_tmp.splice(out_tmp.end(), result.first);;
		}
		assert(result.second.empty());
	}
	if (!out_tmp.empty())
		out.splice(out.end(), recompute_reversed(out_tmp, p));
	return out;
}

list<Hsp> swipe_set(const SequenceSet::ConstIterator begin, const SequenceSet::ConstIterator end, Params& p) {
	const unsigned b = bin(p.v, 0, 0, 0, 0, 0, 0);
	pair<list<Hsp>, TargetVec> result = swipe_bin(b, begin, end, 0, p);
	if (reversed(p.v))
		result.first = recompute_reversed(result.first, p);
	if (b < BINS - 1 && !result.second.empty()) {
		Targets targets;
		targets[b + 1] = std::move(result.second);
		result.first.splice(result.first.end(), swipe(targets, p));
	}
	return result.first;
}

}

DISPATCH_2(std::list<Hsp>, swipe, const Targets&, targets, Params&, params)
DISPATCH_3(std::list<Hsp>, swipe_set, const SequenceSet::ConstIterator, begin, const SequenceSet::ConstIterator, end, Params&, params)
DISPATCH_7(unsigned, bin, HspValues, v, int, query_len, int, score, int, ungapped_score, const int64_t, dp_size, unsigned, score_width, const Loc, mismatch_est)

}}