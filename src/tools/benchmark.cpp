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
#include "../basic/score_matrix.h"
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
#include "../util/simd/transpose16x16.h"
#include "../util/simd/transpose32x32.h"
#include "../dp/scan_diags.h"

using std::vector;
using std::chrono::high_resolution_clock;
using std::chrono::nanoseconds;
using std::chrono::duration_cast;
using std::list;
using namespace DISPATCH_ARCH;

namespace Benchmark { namespace DISPATCH_ARCH {


constexpr int N = 32;

void print(int* m) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j)
			cout << std::setw(5) << *(m++);
		cout << endl;
	}
}

void tran(const int* r1, const int* r2, int* r3, int* r4, int n) {
	for (int i = 0; i < 32; i+=n*2) {
		std::copy(r1, r1 + n, r3);
		std::copy(r2, r2 + n, r3 + n);
		std::copy(r1 + n, r1 + 2 * n, r4);
		std::copy(r2 + n, r2 + 2 * n, r4 + n);
		r1 += 2 * n;
		r2 += 2 * n;
		r3 += 2 * n;
		r4 += 2 * n;
	}
}

void tran2_low(const int* r1, const int* r2, int* r3, int n) {
	for (int i = 0; i < N/2; i += n) {
		std::copy(r1, r1 + n, r3);
		r3 += n;
		std::copy(r2, r2 + n, r3);
		r3 += n;
		r1 += n;
		r2 += n;
	}
}

void tran2_hi(const int* r1, const int* r2, int* r3, int n) {
	r1 += N/2;
	r2 += N/2;
	for (int i = 0; i < N/2; i += n) {
		std::copy(r1, r1 + n, r3);
		r3 += n;
		std::copy(r2, r2 + n, r3);
		r3 += n;
		r1 += n;
		r2 += n;
	}
}

void tran() {
	int m1[32 * 32], m2[32 * 32], m3[32 * 32], m4[32 * 32], m5[32 * 32], m6[32 * 32];
	for (int i = 0; i < 1024; ++i)
		m1[i] = i;
	print(m1);
	for (int i = 0; i < 32; i += 2)
		tran(&m1[32 * i], &m1[32 * (i + 1)], &m2[32 * i], &m2[32 * (i + 1)], 1);
	print(m2);
	for (int i = 0; i < 32; i += 4)
		for (int j = 0; j < 2; ++j)
			tran(&m2[32 * (i+j)], &m2[32 * (i+j + 2)], &m3[32 * (i+j)], &m3[32 * (i+j + 2)], 2);
	print(m3);
	for (int i = 0; i < 32; i += 8)
		for (int j = 0; j < 4; ++j)
			tran(&m3[32 * (i + j)], &m3[32 * (i + j + 4)], &m4[32 * (i + j)], &m4[32 * (i + j + 4)], 4);
	print(m4);
	for (int i = 0; i < 32; i += 16)
		for (int j = 0; j < 8; ++j)
			tran(&m4[32 * (i + j)], &m4[32 * (i + j + 8)], &m5[32 * (i + j)], &m5[32 * (i + j + 8)], 8);
	print(m5);
	for (int j = 0; j < 16; ++j)
		tran(&m5[32 * (j)], &m5[32 * (j + 16)], &m6[32 * (j)], &m6[32 * (j + 16)], 16);
	print(m6);
}

void tran_8x8() {
	int m1[N*N], m2[N * N], m3[N * N], m4[N * N], m5[N * N], m6[N * N];
	for (int i = 0; i < N*N; ++i)
		m1[i] = i;
	print(m1);
	cout << endl;

	for (int i = 0; i < N; i += 2) {
		tran2_low(&m1[N * i], &m1[N * (i + 1)], &m2[N * (i / 2)], 1);
		tran2_hi(&m1[N * i], &m1[N * (i + 1)], &m2[N * (i / 2 + N / 2)], 1);
	}
	print(m2);
	cout << endl;

	for (int i = 0; i < N; i += 2) {
		int block = i / 4;
		int start = block * 4;
		int j = (i % 4) / 2;
		tran2_low(&m2[N * i], &m2[N * (i + 1)], &m3[N * (start + j)], 2);
		tran2_hi(&m2[N * i], &m2[N * (i + 1)], &m3[N * (start + 2 + j)], 2);
	}
	print(m3);
	cout << endl;

	for (int i = 0; i < N; i += 2) {
		int block = i / 2;
		int start = block * 2;
		int j = (i % 2) / 2;
		tran2_low(&m3[N * i], &m3[N * (i + 1)], &m4[N * (start + j)], 4);
		tran2_hi(&m3[N * i], &m3[N * (i + 1)], &m4[N * (start + 1 + j)], 4);
	}
	print(m4);
	cout << endl;
}

void tran_blocks(const int* in, int* out, int block_size, int pack, int *r_in, int *r_out) {
	for (int i = 0; i < N; i += 2) {
		int block = i / block_size;
		int start = block * block_size;
		int j = (i % block_size) / 2;
		tran2_low(&in[N * i], &in[N * (i + 1)], &out[N * (start + j)], pack);
		r_out[start + j] = r_in[i];
		tran2_hi(&in[N * i], &in[N * (i + 1)], &out[N * (start + block_size/2 + j)], pack);
		r_out[start + block_size/2 + j] = r_in[i + 1];
		cout << "t = r" << r_in[i] << ';' << endl;
		cout << "r" << r_in[i] << " = _mm256_unpacklo_epi" << pack * 8 << "(t, r" << r_in[i + 1] << ");" << endl;
		cout << "r" << r_in[i + 1] << " = _mm256_unpackhi_epi" << pack * 8 << "(t, r" << r_in[i + 1] << ");" << endl;
	}
	print(out);
	cout << endl;
}

void tran_16x16() {
	int m1[N * N], m2[N * N], m3[N * N], m4[N * N], m5[N * N], m6[N * N];
	int r1[N], r2[N], r3[N], r4[N], r5[N];
	for (int i = 0; i < N * N; ++i)
		m2[i] = i;
	for (int i = 0; i < N; ++i) {
		cout << "case " << (N-i) << ": r" << i << " = _mm_loadu_si128((const __m256i*)*(data++));" << endl;
		r1[i] = i;
	}
	print(m2);
	cout << endl;
	
	tran_blocks(m2, m3, 16, 1, r1,r2);;
	tran_blocks(m3, m4, 8, 2,r2,r3);
	tran_blocks(m4, m5, 4, 4,r3,r4);
	tran_blocks(m5, m6, 2, 8,r4,r5);

	for (int i = 0; i < 16; ++i)
		cout << "_mm_store_si128(ptr++, r" << r5[i] << ");" << endl;
}

void tran_32x32() {
	int m1[N * N], m2[N * N], m3[N * N], m4[N * N], m5[N * N], m6[N * N];
	int r1[N], r2[N], r3[N], r4[N], r5[N], r6[N];
	for (int i = 0; i < N * N; ++i)
		m1[i] = i;
	for (int i = 0; i < N; ++i)
		cout << "r" << i << ",";
	cout << endl;
	for (int i = 0; i < N; ++i) {
		cout << "case " << (N - i) << ": r" << i << " = _mm256_loadu_si256((const __m256i*)*(data++));" << endl;
		r1[i] = i;
	}
	print(m1);
	cout << endl;

	tran_blocks(m1, m2, 32, 1, r1, r2);
	tran_blocks(m2, m3, 16, 2, r2, r3);
	tran_blocks(m3, m4, 8, 4, r3, r4);
	tran_blocks(m4, m5, 4, 8, r4, r5);
	tran_blocks(m5, m6, 2, 16, r5, r6);

	for (int i = 0; i < N; ++i)
		cout << "_mm256_store_si256(ptr++, r" << r6[i] << ");" << endl;
}

#ifdef __SSE__
void benchmark_transpose2() {
	tran_32x32();
	static const size_t n = 10000000llu;
	static float in[8*8], out[8*8];
	float* m1 = in, * m2 = out;
	//for(int i=0;i<1024;++i)
		//*(char*)

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		//tran(m1, m2);
		std::swap(m1, m2);
		++m1[i & 63];
	}
	cout << "Matrix transpose 32x32 bytes:\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 256) * 1000 << " ps/Letter" << endl;
}
#endif

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
	constexpr size_t CHANNELS = ScoreTraits<score_vector<int8_t>>::CHANNELS;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	const Letter *q = s1.data(), *s = s2.data();
	int score = 0;
	score_vector<int8_t> sv;
	::DISPATCH_ARCH::SIMD::Vector<int8_t> seq(s1.data());

	for (size_t i = 0; i < n; ++i) {		
		sv  = score_vector<int8_t>(i & 15, seq);
		volatile auto x = sv.data_;
	}
	cout << "SSSE3 score shuffle:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * CHANNELS) * 1000 << " ps/Letter" << endl;
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
		::DP::window_ungapped(s1.data(), targets, 16, 64, out);
	}
	cout << "SSE ungapped extend:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 16 * 64) * 1000 << " ps/Cell" << endl;
}
#endif

#ifdef __SSE__
void benchmark_transpose() {
	static const size_t n = 10000000llu;
	static signed char in[256], out[256];
	signed char* v[16];
	for (int i = 0; i < 16; ++i)
		v[i] = &in[i * 16];

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		transpose16x16((const signed char**)v, 16, out);
		in[0] = out[0];
	}
	cout << "Matrix transpose 16x16 bytes:\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 256) * 1000 << " ps/Letter" << endl;

#ifdef __AVX2__
	{
		static signed char in[32 * 32], out[32 * 32];
		signed char* v[32];
		for (int i = 0; i < 32; ++i)
			v[i] = &in[i * 32];

		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		for (size_t i = 0; i < n; ++i) {
			transpose32x32((const signed char**)v, 32, out);
			in[0] = out[0];
		}
		cout << "Matrix transpose 32x32 bytes:\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 32 * 32) * 1000 << " ps/Letter" << endl;
	}
#endif
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
			diagonal_cell = ::swipe_cell_update(diagonal_cell, scores, nullptr, gap_extension, gap_open, horizontal_gap, vertical_gap, best);
		}
		volatile auto x = diagonal_cell.data_;
	}
	cout << "SWIPE cell update (int8_t):\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 16) * 1000 << " ps/Cell" << endl;
#endif
}
#endif

#ifdef __SSE4_1__
void swipe(const sequence &s1, const sequence &s2) {
	constexpr int CHANNELS = ::DISPATCH_ARCH::ScoreTraits<score_vector<int8_t>>::CHANNELS;
	sequence target[CHANNELS];
	std::fill(target, target + CHANNELS, s2);
	static const size_t n = 1000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile list<Hsp> v = ::DP::Swipe::swipe(s1, target, target + CHANNELS, 100);
	}
	cout << "SWIPE (int8_t):\t\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s1.length() * s2.length() * CHANNELS) * 1000 << " ps/Cell" << endl;
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
		volatile auto out = ::DP::BandedSwipe::swipe(s1, target8, target16, Frame(0), &cbs, 0, 0, stat);
	}
	cout << "Banded SWIPE (int16_t, CBS):\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s1.length() * 65 * 8) * 1000 << " ps/Cell" << endl;

	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile auto out = ::DP::BandedSwipe::swipe(s1, target8, target16, Frame(0), nullptr, 0, 0, stat);
	}
	cout << "Banded SWIPE (int16_t):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s1.length() * 65 * 8) * 1000 << " ps/Cell" << endl;
}

#ifdef __SSE4_1__
void diag_scores(const sequence& s1, const sequence& s2) {
	static const size_t n = 100000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	Bias_correction cbs(s1);
	LongScoreProfile p(s1, cbs);
	int scores[128];
	for (size_t i = 0; i < n; ++i) {
		DP::scan_diags128(p, s2, -32, 0, (int)s2.length(), scores);
	}
	cout << "Diagonal scores:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s2.length() * 128) * 1000 << " ps/Cell" << endl;
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

#ifdef __AVX2__
	//benchmark_transpose2();
#endif
#ifdef __SSE4_1__
	swipe(s3, s4);
	diag_scores(s1, s2);
#endif
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
}

}}
