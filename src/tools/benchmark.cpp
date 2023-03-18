/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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

#include <chrono>
#include <utility>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include "../basic/sequence.h"
#include "../stats/score_matrix.h"
#include "../dp/score_vector.h"
#include "../dp/swipe/swipe.h"
#include "../dp/dp.h"
#include "../dp/score_vector_int8.h"
#include "../dp/score_vector_int16.h"
#include "../dp/score_profile.h"
#include "../dp/ungapped.h"
#include "../search/finger_print.h"
#include "../dp/ungapped_simd.h"
#include "../util/simd/vector.h"
#include "../util/simd/transpose.h"
#include "../dp/scan_diags.h"
#include "../stats/cbs.h"
#include "../util/profiler.h"
#include "../dp/swipe/cell_update.h"
#include "../dp/pfscan/pfscan.h"
#include "../dp/pfscan/simd.h"
#include "../dp/swipe/anchored.h"
#include "../dp/swipe/config.h"
#include "../util/simd/dispatch.h"

void benchmark_io();

using std::vector;
using std::endl;
using std::chrono::high_resolution_clock;
using std::chrono::nanoseconds;
using std::chrono::duration_cast;
using std::list;
using std::array;
using namespace DISPATCH_ARCH;

namespace Benchmark { namespace DISPATCH_ARCH {

#if defined(__SSE4_1__) && defined(EXTRA)
void swipe_cell_update();
#endif

#ifdef __SSE4_1__
void benchmark_hamming(const Sequence& s1, const Sequence& s2) {
	static const size_t n = 100000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	Byte_finger_print_48 f1(s1.data()), f2 (s2.data());
	for (size_t i = 0; i < n; ++i) {
		f1.r1 = _mm_xor_si128(f1.r1, f1.r2);
		volatile unsigned y = f1.match(f2);
	}
	message_stream << "SSE hamming distance:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 48) * 1000 << " ps/Cell" << endl;
}
#endif

void benchmark_ungapped(const Sequence& s1, const Sequence& s2)
{
	static const size_t n = 10000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	const Letter *q = s1.data(), *s = s2.data();

	for (size_t i = 0; i < n; ++i) {

		volatile int score = ungapped_window(q, s, 64);

	}
	
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	std::chrono::nanoseconds time_span = duration_cast<std::chrono::nanoseconds>(t2 - t1);

	message_stream << "Scalar ungapped extension:\t" << (double)time_span.count() / (n*64) * 1000 << " ps/Cell" << endl;
}

#if defined(__SSSE3__) && defined(__SSE4_1__)
void benchmark_ssse3_shuffle(const Sequence&s1, const Sequence&s2)
{
	static const size_t n = 100000000llu;
	constexpr size_t CHANNELS = ScoreTraits<ScoreVector<int8_t, SCHAR_MIN>>::CHANNELS;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	const Letter *q = s1.data(), *s = s2.data();
	int score = 0;
	ScoreVector<int8_t, SCHAR_MIN> sv;
	::DISPATCH_ARCH::SIMD::Vector<int8_t> seq(s1.data());

	for (size_t i = 0; i < n; ++i) {		
		sv  = ScoreVector<int8_t, SCHAR_MIN>(i & 15, seq);
		volatile auto x = sv.data_;
	}
	message_stream << "SSSE3 score shuffle:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * CHANNELS) * 1000 << " ps/Letter" << endl;
}
#endif

#ifdef __SSE4_1__
void benchmark_ungapped_sse(const Sequence&s1, const Sequence&s2) {
	static const size_t n = 1000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	const Letter* targets[16];
	//int out[16];
	for (int i = 0; i < 16; ++i)
		targets[i] = s2.data();

	for (size_t i = 0; i < n; ++i) {
		//::DP::ARCH_SSE4_1::window_ungapped(s1.data(), targets, 16, 64, out);
	}
	message_stream << "SSE ungapped extend:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 16 * 64) * 1000 << " ps/Cell" << endl;

#ifdef __AVX2__
	{
		high_resolution_clock::time_point t1 = high_resolution_clock::now();

		const Letter* targets[32];
		// int out[32];
		for (int i = 0; i < 32; ++i)
			targets[i] = s2.data();

		for (size_t i = 0; i < n; ++i) {
			//::DP::ARCH_AVX2::window_ungapped(s1.data(), targets, 32, 64, out);
		}
		message_stream << "AVX2 ungapped extend:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 32 * 64) * 1000 << " ps/Cell" << endl;
	}
#endif
}
#endif

#ifdef __SSE2__
void benchmark_transpose() {
	static const size_t n = 10000000llu;
	static signed char in[256], out[256];
	signed char* v[16];
	for (int i = 0; i < 16; ++i)
		v[i] = &in[i * 16];

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		transpose((const signed char**)v, 16, out, __m128i());
		in[0] = out[0];
	}
	message_stream << "Matrix transpose 16x16 bytes:\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 256) * 1000 << " ps/Letter" << endl;

#if ARCH_ID == 2
	{
		static signed char in[32 * 32], out[32 * 32];
		signed char* v[32];
		for (int i = 0; i < 32; ++i)
			v[i] = &in[i * 32];

		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		for (size_t i = 0; i < n; ++i) {
			transpose((const signed char**)v, 32, out, __m256i());
			in[0] = out[0];
		}
		message_stream << "Matrix transpose 32x32 bytes:\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 32 * 32) * 1000 << " ps/Letter" << endl;
	}
#endif
}
#endif

#ifdef __SSE4_1__

void mt_swipe(const Sequence& s1, const Sequence& s2) {
	//constexpr int CHANNELS = 16;
	constexpr int CHANNELS = ::DISPATCH_ARCH::ScoreTraits<ScoreVector<int8_t, SCHAR_MIN>>::CHANNELS;
	static const size_t n = 100000llu;
	DP::Targets targets;
	for (size_t i = 0; i < CHANNELS; ++i)
		targets[0].emplace_back(s2, s2.length(), 0, 0, Interval(), 0, 0, 0);
	Bias_correction cbs(s1);
	Statistics stat;
	Sequence query = s1;
	query.len_ = std::min(query.len_, (Loc)255);
	auto dp_size = (n * query.length() * s2.length() * CHANNELS);
	DP::Params params{
		query, "", Frame(0), query.length(), cbs.int8.data(), DP::Flags::FULL_MATRIX, HspValues(), stat, nullptr
	};

	auto f = [&]() {
		for (size_t i = 0; i < n; ++i) {
			//volatile list<Hsp> v = ::DP::BandedSwipe::ARCH_SSE4_1::swipe(targets, params);
			volatile list<Hsp> v = ::DP::BandedSwipe::swipe(targets, params);
		}
	};
	using std::thread;
	vector<thread> th;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (int i = 0; i < config.threads_; ++i)
		th.emplace_back(f);
	for (auto& t : th)
		t.join();
	message_stream << "MT_SWIPE (int8_t):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / dp_size * 1000 << " ps/Cell" << endl;
}

void swipe(const Sequence&s1, const Sequence&s2) {
	constexpr int CHANNELS = ::DISPATCH_ARCH::ScoreTraits<ScoreVector<int8_t, SCHAR_MIN>>::CHANNELS;
	static const size_t n = 1000llu;
	DP::Targets targets;
	for (size_t i = 0; i < 32; ++i)
		targets[0].emplace_back(s2, s2.length(), 0, 0, Interval(), 0, 0, 0);
	Bias_correction cbs(s1);
	Statistics stat;
	Sequence query = s1;
	query.len_ = std::min(query.len_, (Loc)255);
	auto dp_size = (n * query.length() * s2.length() * CHANNELS);
	config.comp_based_stats = 4;
	Stats::TargetMatrix matrix(Stats::composition(s1), s1.length(), s2);
	DP::Params params{
		query, "", Frame(0), query.length(), cbs.int8.data(), DP::Flags::FULL_MATRIX, HspValues(), stat, nullptr
	};

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile list<Hsp> v = ::DP::BandedSwipe::swipe(targets, params);
	}
	message_stream << "SWIPE (int8_t):\t\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / dp_size * 1000 << " ps/Cell" << endl;

	t1 = high_resolution_clock::now();
	targets[1] = targets[0];
	targets[0].clear();
	for (size_t i = 0; i < n; ++i) {
		volatile list<Hsp> v = ::DP::BandedSwipe::swipe(targets, params);
	}
	message_stream << "SWIPE (int16_t):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / dp_size * 1000 << " ps/Cell" << endl;

	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile list<Hsp> v = ::DP::BandedSwipe::swipe(targets, params);
	}
	message_stream << "SWIPE (int8_t, Stats):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / dp_size * 1000 << " ps/Cell" << endl;

	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < 32; ++i)
		targets[0][i].matrix = &matrix;
	for (size_t i = 0; i < n; ++i) {
		volatile list<Hsp> v = ::DP::BandedSwipe::swipe(targets, params);
	}
	message_stream << "SWIPE (int8_t, MatrixAdjust):\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / dp_size * 1000 << " ps/Cell" << endl;

	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile list<Hsp> v = ::DP::BandedSwipe::swipe(targets, params);
	}
	message_stream << "SWIPE (int8_t, CBS):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / dp_size * 1000 << " ps/Cell" << endl;

	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile list<Hsp> v = ::DP::BandedSwipe::swipe(targets, params);
	}
	message_stream << "SWIPE (int8_t, TB):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / dp_size * 1000 << " ps/Cell" << endl;
}
#endif

void banded_swipe(const Sequence &s1, const Sequence &s2) {
	DP::Targets targets;
	//config.traceback_mode = TracebackMode::SCORE_BUFFER;
	for (size_t i = 0; i < 8; ++i)
		targets[1].emplace_back(s2, s2.length(), -32, 32, Interval(), 0, 0, 0);
	static const size_t n = 10000llu;
	//static const size_t n = 1llu;
	Statistics stat;
	Bias_correction cbs(s1);
	DP::Params params{
		s1, "", Frame(0), s1.length(), cbs.int8.data(), DP::Flags::NONE, HspValues(), stat, nullptr
	};
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile auto out = ::DP::BandedSwipe::swipe(targets, params);
	}
	message_stream << "Banded SWIPE (int16_t, CBS):\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s1.length() * 65 * 16) * 1000 << " ps/Cell" << endl;
	
	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile auto out = ::DP::BandedSwipe::swipe(targets, params);
	}
	message_stream << "Banded SWIPE (int16_t):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s1.length() * 65 * 16) * 1000 << " ps/Cell" << endl;

	params.v = HspValues::TRANSCRIPT;
	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile auto out = ::DP::BandedSwipe::swipe(targets, params);
	}
	message_stream << "Banded SWIPE (int16_t, CBS, TB):" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s1.length() * 65 * 16) * 1000 << " ps/Cell" << endl;
}

#if ARCH_ID == 2

void anchored_swipe(const Sequence& s1, const Sequence& s2) {
	static const size_t n = 10000llu;
	const auto s1_ = s1.subseq(0, 128);
	const auto s2_ = s2.subseq(0, 128);
	LongScoreProfile<int16_t> prof = DP::make_profile16(s1_, nullptr, 0);
	LongScoreProfile<int8_t> prof8 = DP::make_profile8(s1_, nullptr, 0);
	auto v = prof.pointers(0);
	auto v8 = prof8.pointers(0);
	vector<DP::AnchoredSwipe::Target<int8_t>> targets;
	vector<LongScoreProfile<int8_t>> profiles(32, prof8);
	vector<vector<const int8_t*>> pointers;
	Statistics stats;
	DP::AnchoredSwipe::Options options{ v.data(), v.data() };
	for (int i = 0; i < 32; ++i) {
		pointers.push_back(profiles[i].pointers(0));
		//targets.push_back(DP::AnchoredSwipe::Target<int8_t> {s2_, -32, 32, { pointers.back().data(), nullptr}, s1_.length() });
		//targets.push_back(DP::AnchoredSwipe::Target<int8_t>(s2_, -32, 32, pointers[0].data(), s1_.length(), 0, false));
		targets.push_back(DP::AnchoredSwipe::Target<int8_t>(s2_, -32, 32, 0, s1_.length(), 0, false));
	}
	const int cols = round_up(s2_.length(), DP::AnchoredSwipe::ARCH_AVX2::L);

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		DP::AnchoredSwipe::ARCH_AVX2::smith_waterman<ScoreVector<int8_t, 0>>(targets.data(), 32, options);
		volatile auto x = targets[0].score;
	}
	message_stream << "Anchored Swipe (int8_t):\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * cols * 64 * 32) * 1000 << " ps/Cell" << endl;

	vector<DP::AnchoredSwipe::Target<int16_t>> targets16;
	for (int i = 0; i < 16; ++i) {
		//targets16.push_back(DP::AnchoredSwipe::Target<int16_t>(s2_, -32, 32, v.data(), s1_.length(), 0, false));
		targets16.push_back(DP::AnchoredSwipe::Target<int16_t>(s2_, -32, 32, 0, s1_.length(), 0, false));
	}
	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		DP::AnchoredSwipe::ARCH_AVX2::smith_waterman<ScoreVector<int16_t, 0>>(targets16.data(), 16, options);
		volatile auto x = targets[0].score;
	}
	message_stream << "Anchored Swipe (int16_t):\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * cols * 64 * 16) * 1000 << " ps/Cell" << endl;

	DP::Targets dp_targets;
	Anchor a(DiagonalSegment(0, 0, 0, 0), 0, 0, 0, 0, 0);
	for (int i = 0; i < 16; ++i)
		dp_targets[0].emplace_back(s2_, s2_.length(), -32, 32, Interval(), 0, 0, s1_.length(), nullptr, DpTarget::CarryOver(), a);
	DP::AnchoredSwipe::Config cfg{ s1_, nullptr, 0, stats, nullptr };

	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		DP::BandedSwipe::anchored_swipe(dp_targets, cfg);
		volatile auto x = targets[0].score;
	}
	message_stream << "Anchored Swipe2 (int16_t):\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 128 * 64 * 16) * 1000 << " ps/Cell" << endl;
}

//#endif
#endif

#ifdef __SSE2__

void prefix_scan(const Sequence& s1, const Sequence& s2) {
	static const size_t n = 100000llu;
	const auto s1_ = s1.subseq(0, 150);
	const auto s2_ = s2.subseq(0, 150);
	Statistics stat;
	Bias_correction cbs(s1_);
	LongScoreProfile<int16_t> prof = DP::make_profile16(s1_, nullptr, 0);
	LongScoreProfile<int8_t> prof8 = DP::make_profile8(s1_, nullptr, 0);
	auto v = prof.pointers(0);
	auto v8 = prof8.pointers(0);
	DP::PrefixScan::Config cfg{ s1_,s2_, "", "", -32, 32, Interval(), v.data(), nullptr, v8.data(), nullptr, stat, 1000, 60, 0 };
	
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		//volatile auto out = DP::PrefixScan::align16(cfg);
	}
	message_stream << "Prefix Scan (int16_t):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (s2_.length() * 64 * n) * 1000 << " ps/Cell" << endl;
	//statistics += stat;
	//statistics.print();
	//statistics.reset();

	cfg.score_bias = 0;
	stat.reset();
	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		//volatile auto out = DP::PrefixScan::align8(cfg);
	}
	message_stream << "Prefix Scan (int8_t):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (s2_.length() * 64 * n) * 1000 << " ps/Cell" << endl;
	statistics += stat;
	//statistics.print();
}

#endif

#ifdef __SSE4_1__
void diag_scores(const Sequence& s1, const Sequence& s2) {
	static const size_t n = 100000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	Bias_correction cbs(s1);
	LongScoreProfile<int8_t> p = DP::make_profile8(s1, cbs.int8.data(), 0);
	int scores[128];
	for (size_t i = 0; i < n; ++i) {
		DP::scan_diags128(p, s2, -32, 0, (int)s2.length(), scores);
		volatile int x = scores[i & 128];
	}
	message_stream << "Diagonal scores:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s2.length() * 128) * 1000 << " ps/Cell" << endl;
}
#endif

void evalue() {
	static const size_t n = 1000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	volatile double x = 0.0;
	for (size_t i = 0; i < n; ++i) {
		x += score_matrix.evalue_norm((int)i, 300);
	}
	message_stream << "Evalue:\t\t\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n) << " ns" << endl;

	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		x += score_matrix.evalue(300, 300, 300);
	}
	message_stream << "Evalue (ALP):\t\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n) << " ns" << endl;
}

void matrix_adjust(const Sequence& s1, const Sequence& s2) {
	static const size_t n = 10000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	vector<double> mat_final(TRUE_AA * TRUE_AA);
	int iteration_count;
	const double* joint_probs = (const double*)(Stats::blosum62.joint_probs);
	auto row_probs = Stats::composition(s1), col_probs = Stats::composition(s2);
	config.cbs_err_tolerance = 0.0001;

	for (size_t i = 0; i < n; ++i) {
		Stats::Blast_OptimizeTargetFrequencies(mat_final.data(),
			TRUE_AA,
			&iteration_count,
			joint_probs,
			row_probs.data(), col_probs.data(),
			true,
			0.44,
			config.cbs_err_tolerance,
			config.cbs_it_limit);
	}

	message_stream << "Matrix adjust:\t\t\t" << (double)duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t1).count() / (n) << " ms" << endl;

	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		Stats::OptimizeTargetFrequencies(mat_final.data(), joint_probs, row_probs.data(), col_probs.data(), 0.44, config.cbs_err_tolerance, config.cbs_it_limit);
	}

	message_stream << "Matrix adjust (vectorized):\t" << (double)duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t1).count() / (n) << " micros" << endl;

	//Profiler::print(n);
}

void benchmark() {
	if (config.type == "swipe") {
#if defined(__SSE4_1__) && defined(EXTRA)
		swipe_cell_update();
#endif
		return;
	}if (!config.type.empty()) {
		benchmark_io();
		return;
	}

	vector<Letter> s1, s2, s3, s4;
		
	s1 = Sequence::from_string("mpeeeysefkelilqkelhvvyalshvcgqdrtllasillriflhekleslllctlndreismedeattlfrattlastlmeqymkatatqfvhhalkdsilkimeskqscelspskleknedvntnlthllnilselvekifmaseilpptlryiygclqksvqhkwptnttmrtrvvsgfvflrlicpailnprmfniisdspspiaartlilvaksvqnlanlvefgakepymegvnpfiksnkhrmimfldelgnvpelpdttehsrtdlsrdlaalheicvahsdelrtlsnergaqqhvlkkllaitellqqkqnqyt"); // d1wera_
	s2 = Sequence::from_string("erlvelvtmmgdqgelpiamalanvvpcsqwdelarvlvtlfdsrhllyqllwnmfskeveladsmqtlfrgnslaskimtfcfkvygatylqklldpllrivitssdwqhvsfevdptrlepsesleenqrnllqmtekffhaiissssefppqlrsvchclyqvvsqrfpqnsigavgsamflrfinpaivspyeagildkkpppiierglklmskilqsianhvlftkeehmrpfndfvksnfdaarrffldiasdcptsdavnhslsfisdgnvlalhrllwnnqekigqylssnrdhkavgrrpfdkmatllaylgppe"); // d1nf1a_
	s3 = Sequence::from_string("ttfgrcavksnqagggtrshdwwpcqlrldvlrqfqpsqnplggdfdyaeafqsldyeavkkdiaalmtesqdwwpadfgnygglfvrmawhsagtyramdgrggggmgqqrfaplnswpdnqnldkarrliwpikqkygnkiswadlmlltgnvalenmgfktlgfgggradtwqsdeavywgaettfvpqgndvrynnsvdinaradklekplaathmgliyvnpegpngtpdpaasakdireafgrmgmndtetvaliagghafgkthgavkgsnigpapeaadlgmqglgwhnsvgdgngpnqmtsgleviwtktptkwsngyleslinnnwtlvespagahqweavngtvdypdpfdktkfrkatmltsdlalindpeylkisqrwlehpeeladafakawfkllhrdlgpttrylgpevp"); // d3ut2a1
	s4 = Sequence::from_string("lvhvasvekgrsyedfqkvynaialklreddeydnyigygpvlvrlawhisgtwdkhdntggsyggtyrfkkefndpsnaglqngfkflepihkefpwissgdlfslggvtavqemqgpkipwrcgrvdtpedttpdngrlpdadkdagyvrtffqrlnmndrevvalmgahalgkthlknsgyegpggaannvftnefylnllnedwklekndanneqwdsksgymmlptdysliqdpkylsivkeyandqdkffkdfskafekllengitfpkdapspfifktleeqgl"); // d2euta_

	Sequence ss1 = Sequence(s1).subseq(34, (Loc)s1.size());
	Sequence ss2 = Sequence(s2).subseq(33, (Loc)s2.size());

#ifdef __SSE4_1__
	//mt_swipe(s3, s4);
#endif
#if ARCH_ID == 2
//#ifdef __SSE4_1__
	anchored_swipe(s1, s2);
	//minimal_sw(s1, s2);
//#endif
#endif
#ifdef __SSE2__
	prefix_scan(s1, s2);
#endif
#ifdef __SSE4_1__
	swipe(s3, s4);
	diag_scores(s1, s2);
#endif
#ifdef __SSE2__
	banded_swipe(s1, s2);
#endif
	evalue();
	matrix_adjust(s1, s2);
#ifdef __SSE4_1__
	benchmark_hamming(s1, s2);
#endif
	benchmark_ungapped(ss1, ss2);
#if defined(__SSSE3__) && defined(__SSE4_1__)
	benchmark_ssse3_shuffle(s1, s2);
#endif
#ifdef __SSE4_1__
	benchmark_ungapped_sse(ss1, ss2);
#endif
#ifdef __SSE2__
	benchmark_transpose();
#endif
}

}

DISPATCH_0V(benchmark)

}
