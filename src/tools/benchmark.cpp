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
#include "../basic/sequence.h"
#include "../basic/score_matrix.h"
#include "../dp/score_vector.h"
#include "../util/simd/transpose.h"
#include "../dp/swipe/swipe.h"
#include "../dp/dp.h"
#include "../dp/score_vector_int8.h"
#include "../dp/score_profile.h"
#include "../dp/ungapped.h"
#include "../search/finger_print.h"
#include "../dp/ungapped_simd.h"

using std::vector;
using std::chrono::high_resolution_clock;
using std::chrono::nanoseconds;
using std::chrono::duration_cast;
using std::list;

namespace Benchmark { namespace DISPATCH_ARCH {

static int xdrop_window(const Letter *query, const Letter *subject) {
	static const int window = 64;
	int score(0), st(0), n = 0;
	const Letter *q(query), *s(subject);

	st = score;
	while (n < window)
	{
		st += score_matrix(*q, *s);
		score = std::max(score, st);
		++q;
		++s;
		++n;

		st += score_matrix(*q, *s);
		score = std::max(score, st);
		++q;
		++s;
		++n;

		st += score_matrix(*q, *s);
		score = std::max(score, st);
		++q;
		++s;
		++n;

		st += score_matrix(*q, *s);
		score = std::max(score, st);
		++q;
		++s;
		++n;
	}
	return score;
}

#ifdef __SSE4_1__
void benchmark_hamming(const sequence &s1, const sequence &s2) {
	static const size_t n = 100000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	Byte_finger_print_48 f1(s1.data()), f2 (s2.data());
	for (size_t i = 0; i < n; ++i) {
		f1.r1 = _mm_xor_si128(f1.r1, f1.r2);
		volatile unsigned y = f1.match(f2);
	}
	cout << "SSE hamming distance:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 48) * 1000 << " ps/Cell" << endl;
}
#endif

void benchmark_ungapped(const sequence &s1, const sequence &s2)
{
	static const size_t n = 10000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	const Letter *q = s1.data(), *s = s2.data();

	for (size_t i = 0; i < n; ++i) {

		volatile int score = xdrop_window(q, s);

	}
	
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	std::chrono::nanoseconds time_span = duration_cast<std::chrono::nanoseconds>(t2 - t1);

	cout << "Scalar ungapped extension:\t" << (double)time_span.count() / (n*64) * 1000 << " ps/Cell" << endl;
}

#ifdef __SSSE3__
void benchmark_ssse3_shuffle(const sequence &s1, const sequence &s2)
{
	static const size_t n = 100000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	const Letter *q = s1.data(), *s = s2.data();
	int score = 0;
	score_vector<uint8_t> sv;
	__m128i seq = _mm_loadu_si128((const __m128i*)s1.data());

	for (size_t i = 0; i < n; ++i) {		
		sv.set_ssse3(i & 15, seq);
		volatile __m128i x = sv.data_;
	}
	cout << "SSSE3 score shuffle:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 16) * 1000 << " ps/Letter" << endl;
}
#endif

#ifdef __SSE4_1__
void benchmark_ungapped_sse(const sequence &s1, const sequence &s2) {
	static const size_t n = 1000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	const Letter* targets[16];
	int out[16];
	for (int i = 0; i < 16; ++i)
		targets[i] = s2.data();

	for (size_t i = 0; i < n; ++i) {
		DP::window_ungapped(s1.data(), targets, 16, 64, out);
	}
	cout << "SSE ungapped extend:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 16 * 64) * 1000 << " ps/Cell" << endl;
}
#endif

#ifdef __SSE__
void benchmark_transpose() {
	static const size_t n = 10000000llu;
	static char in[256], out[256];

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		transpose(in, out, 0);
		in[0] = out[0];
	}
	cout << "Matrix transpose 16x16 bytes:\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 256) * 1000 << " ps/Letter" << endl;
}
#endif

#ifdef __SSE2__
void swipe_cell_update() {
	static const size_t n = 1000000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	{
		score_vector<uint8_t> diagonal_cell, scores, gap_extension, gap_open, horizontal_gap, vertical_gap, best, vbias;
		for (size_t i = 0; i < n; ++i) {
			diagonal_cell = cell_update(diagonal_cell, scores, gap_extension, gap_open, horizontal_gap, vertical_gap, best, vbias);
		}
		volatile __m128i x = diagonal_cell.data_;
	}
	cout << "SWIPE cell update (uint8_t):\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 16) * 1000 << " ps/Cell" << endl;

#ifdef __SSE4_1__
	t1 = high_resolution_clock::now();
	{
		score_vector<int8_t> diagonal_cell, scores, gap_extension, gap_open, horizontal_gap, vertical_gap, best;
		for (size_t i = 0; i < n; ++i) {
			diagonal_cell = swipe_cell_update(diagonal_cell, scores, nullptr, gap_extension, gap_open, horizontal_gap, vertical_gap, best);
		}
		volatile auto x = diagonal_cell.data_;
	}
	cout << "SWIPE cell update (int8_t):\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 16) * 1000 << " ps/Cell" << endl;
#endif
}
#endif

#ifdef __SSE4_1__
void swipe(const sequence &s1, const sequence &s2) {
	sequence target[16];
	std::fill(target, target + 16, s2);
	static const size_t n = 10000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile list<Hsp> v = DP::Swipe::swipe(s1, target, target + 16, 100);
	}
	cout << "SWIPE (int8_t):\t\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s1.length() * s2.length() * 16) * 1000 << " ps/Cell" << endl;
}
#endif

void banded_swipe(const sequence &s1, const sequence &s2) {
	vector<DpTarget> target8, target16;
	for (size_t i = 0; i < 8; ++i)
		target16.emplace_back(s2, -32, 32, 0, 0);
	static const size_t n = 10000llu;
	//static const size_t n = 1llu;
	Statistics stat;
	Bias_correction cbs(s1);
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile auto out = DP::BandedSwipe::swipe(s1, target8, target16, Frame(0), &cbs, 0, 0, stat);
	}
	cout << "Banded SWIPE (int16_t, CBS):\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s1.length() * 65 * 8) * 1000 << " ps/Cell" << endl;

	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile auto out = DP::BandedSwipe::swipe(s1, target8, target16, Frame(0), nullptr, 0, 0, stat);
	}
	cout << "Banded SWIPE (int16_t):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s1.length() * 65 * 8) * 1000 << " ps/Cell" << endl;
}

#ifdef __SSE4_1__
void diag_scores(const sequence &s1, const sequence &s2) {
	static const size_t n = 100000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	Bias_correction cbs(s1);
	LongScoreProfile p(s1, cbs);
	for (size_t i = 0; i < n; ++i) {
		//scan_cols(p, s2, 0, 0, (int)s2.length());
	}
	cout << "Diagonal scores:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s2.length() * 64) * 1000 << " ps/Cell" << endl;
}
#endif

#ifdef __SSE4_1__
void diag_scores2(const sequence& s1, const sequence& s2) {
	static const size_t n = 100000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	Bias_correction cbs(s1);
	LongScoreProfile p(s1, cbs);
	int scores[64];
	for (size_t i = 0; i < n; ++i) {
		DP::scan_diags128(p, s2, -32, 0, (int)s2.length(), &scores[0]);
	}
	cout << "Diagonal scores (2):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s2.length() * 128) * 1000 << " ps/Cell" << endl;
}
#endif

void benchmark() {
	vector<Letter> s1, s2, s3, s4;
		
	s1 = sequence::from_string("mpeeeysefkelilqkelhvvyalshvcgqdrtllasillriflhekleslllctlndreismedeattlfrattlastlmeqymkatatqfvhhalkdsilkimeskqscelspskleknedvntnlthllnilselvekifmaseilpptlryiygclqksvqhkwptnttmrtrvvsgfvflrlicpailnprmfniisdspspiaartlilvaksvqnlanlvefgakepymegvnpfiksnkhrmimfldelgnvpelpdttehsrtdlsrdlaalheicvahsdelrtlsnergaqqhvlkkllaitellqqkqnqyt"); // d1wera_
	s2 = sequence::from_string("erlvelvtmmgdqgelpiamalanvvpcsqwdelarvlvtlfdsrhllyqllwnmfskeveladsmqtlfrgnslaskimtfcfkvygatylqklldpllrivitssdwqhvsfevdptrlepsesleenqrnllqmtekffhaiissssefppqlrsvchclyqvvsqrfpqnsigavgsamflrfinpaivspyeagildkkpppiierglklmskilqsianhvlftkeehmrpfndfvksnfdaarrffldiasdcptsdavnhslsfisdgnvlalhrllwnnqekigqylssnrdhkavgrrpfdkmatllaylgppe"); // d1nf1a_
	s3 = sequence::from_string("ttfgrcavksnqagggtrshdwwpcqlrldvlrqfqpsqnplggdfdyaeafqsldyeavkkdiaalmtesqdwwpadfgnygglfvrmawhsagtyramdgrggggmgqqrfaplnswpdnqnldkarrliwpikqkygnkiswadlmlltgnvalenmgfktlgfgggradtwqsdeavywgaettfvpqgndvrynnsvdinaradklekplaathmgliyvnpegpngtpdpaasakdireafgrmgmndtetvaliagghafgkthgavkgsnigpapeaadlgmqglgwhnsvgdgngpnqmtsgleviwtktptkwsngyleslinnnwtlvespagahqweavngtvdypdpfdktkfrkatmltsdlalindpeylkisqrwlehpeeladafakawfkllhrdlgpttrylgpevp"); // d3ut2a1
	s4 = sequence::from_string("lvhvasvekgrsyedfqkvynaialklreddeydnyigygpvlvrlawhisgtwdkhdntggsyggtyrfkkefndpsnaglqngfkflepihkefpwissgdlfslggvtavqemqgpkipwrcgrvdtpedttpdngrlpdadkdagyvrtffqrlnmndrevvalmgahalgkthlknsgyegpggaannvftnefylnllnedwklekndanneqwdsksgymmlptdysliqdpkylsivkeyandqdkffkdfskafekllengitfpkdapspfifktleeqgl"); // d2euta_

	sequence ss1 = sequence(s1).subseq(34, s1.size());
	sequence ss2 = sequence(s2).subseq(33, s2.size());

#ifdef __SSE4_1__
	benchmark_hamming(s1, s2);
#endif
	benchmark_ungapped(ss1, ss2);
#ifdef __SSSE3__
	benchmark_ssse3_shuffle(s1, s2);
#endif
#ifdef __SSE4_1__
	benchmark_ungapped_sse(ss1, ss2);
#endif
#ifdef __SSE__
	benchmark_transpose();
#endif
#ifdef __SSE2__
	banded_swipe(s1, s2);
	swipe_cell_update();
#endif
#ifdef __SSE4_1__
	swipe(s3, s4);
	diag_scores(s1, s2);
	diag_scores2(s1, s2);
#endif
}

}}
