#include <chrono>
#include <utility>
#include "../basic/sequence.h"
#include "../basic/score_matrix.h"
#include "../dp/score_vector.h"
#include "../util/simd/transpose.h"
#include "../dp/swipe/swipe.h"
#include "../dp/dp.h"

using std::vector;
using std::chrono::high_resolution_clock;
using std::chrono::nanoseconds;
using std::chrono::duration_cast;
using std::list;

namespace Benchmark { namespace DISPATCH_ARCH {

void scan_cols(const Long_score_profile &qp, sequence s, int i, int j, int j_end)
{
#ifdef __SSE2__
	typedef score_vector<int8_t> Sv;
	const int qlen = (int)qp.length();

	int j2 = std::max(-(i - j + 15), j),
		i3 = j2 + i - j,
		j2_end = std::min(qlen - (i - j), j_end);
	Sv v, max, v1, max1, v2, max2, v3, max3;
	for (; j2 < j2_end; ++j2, ++i3) {
		const int8_t *q = (int8_t*)qp.get(s[j2], i3);
		v = v + score_vector<int8_t>((int8_t*)q);
		max.max(v);
		q += 16;
		v1 = v1 + score_vector<int8_t>((int8_t*)q);
		max1.max(v1);
		q += 16;
		v2 = v2 + score_vector<int8_t>((int8_t*)q);
		max2.max(v2);
		q += 16;
		v3 = v3 + score_vector<int8_t>((int8_t*)q);
		max3.max(v3);
	}
	volatile __m128i x = max.data_, y = max1.data_, z = max2.data_, w = max3.data_;
#endif
}

int xdrop_window2(const Letter *query, const Letter *subject)
{
	static const int window = 40;
	int score(0), st(0), n = 0, i = 0;

	const Letter *q(query), *s(subject);

	st = score;
	while (n < window)
	{
		st += score_matrix(*q, *s);
		if (st > score) {
			score = st;
			i = n;
		}
		++q;
		++s;
		++n;

		st += score_matrix(*q, *s);
		if (st > score) {
			score = st;
			i = n;
		}
		++q;
		++s;
		++n;

		st += score_matrix(*q, *s);
		if (st > score) {
			score = st;
			i = n;
		}
		++q;
		++s;
		++n;
	}
	return st * i;
}

void benchmark_ungapped(const sequence &s1, const sequence &s2)
{
	static const size_t n = 10000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	const Letter *q = s1.data(), *s = s2.data();

	for (size_t i = 0; i < n; ++i) {

		volatile int score = xdrop_window2(q, s);

	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	std::chrono::nanoseconds time_span = duration_cast<std::chrono::nanoseconds>(t2 - t1);

	cout << "Scalar ungapped extension:\t" << (double)time_span.count() / (n*40) * 1000 << " ps/Cell" << endl;
}

#ifdef __SSSE3__
void benchmark_ungapped_sse(const sequence &s1, const sequence &s2)
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
	cout << "SSE score shuffle:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 16) * 1000 << " ps/Letter" << endl;
}
#endif

void benchmark_transpose() {
	static const size_t n = 100000000llu;
	static char in[256], out[256];

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		transpose(in, out, 0);
		in[0] = out[0];
	}
	cout << "Matrix transpose 16x16 bytes:\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 256) * 1000 << " ps/Letter" << endl;
}

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
			diagonal_cell = cell_update_sv(diagonal_cell, scores, gap_extension, gap_open, horizontal_gap, vertical_gap, best);
		}
		volatile __m128i x = diagonal_cell.data_;
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
	vector<DpTarget> target;
	for (size_t i = 0; i < 8; ++i)
		target.emplace_back(s2, -32, 32);
	static const size_t n = 10000llu;
	Bias_correction cbs(s1);
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile auto out = DP::BandedSwipe::swipe(s1, target.begin(), target.end(), Frame(0), &cbs, 0, 0);
	}
	cout << "Banded SWIPE (CBS):\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s1.length() * 65 * 8) * 1000 << " ps/Cell" << endl;

	t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		volatile auto out = DP::BandedSwipe::swipe(s1, target.begin(), target.end(), Frame(0), nullptr, 0, 0);
	}
	cout << "Banded SWIPE:\t\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s1.length() * 65 * 8) * 1000 << " ps/Cell" << endl;
}

#ifdef __SSE2__
void diag_scores(const sequence &s1, const sequence &s2) {
	static const size_t n = 100000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	Long_score_profile p(s1);
	for (size_t i = 0; i < n; ++i) {
		scan_cols(p, s2, 0, 0, (int)s2.length());
	}
	cout << "Diagonal scores:\t\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * s2.length() * 64) * 1000 << " ps/Cell" << endl;
}
#endif

void benchmark() {
	vector<Letter> s1, s2, s3, s4;
	
	s1 = sequence::from_string("mpeeeysefkelilqkelhvvyalshvcgqdrtllasillriflhekleslllctlndreismedeattlfrattlastlmeqymkatatqfvhhalkdsilkimeskqscelspskleknedvntnlthllnilselvekifmaseilpptlryiygclqksvqhkwptnttmrtrvvsgfvflrlicpailnprmfniisdspspiaartlilvaksvqnlanlvefgakepymegvnpfiksnkhrmimfldelgnvpelpdttehsrtdlsrdlaalheicvahsdelrtlsnergaqqhvlkkllaitellqqkqnqyt");
	s2 = sequence::from_string("erlvelvtmmgdqgelpiamalanvvpcsqwdelarvlvtlfdsrhllyqllwnmfskeveladsmqtlfrgnslaskimtfcfkvygatylqklldpllrivitssdwqhvsfevdptrlepsesleenqrnllqmtekffhaiissssefppqlrsvchclyqvvsqrfpqnsigavgsamflrfinpaivspyeagildkkpppiierglklmskilqsianhvlftkeehmrpfndfvksnfdaarrffldiasdcptsdavnhslsfisdgnvlalhrllwnnqekigqylssnrdhkavgrrpfdkmatllaylgppe");
	s3 = sequence::from_string("ttfgrcavksnqagggtrshdwwpcqlrldvlrqfqpsqnplggdfdyaeafqsldyeavkkdiaalmtesqdwwpadfgnygglfvrmawhsagtyramdgrggggmgqqrfaplnswpdnqnldkarrliwpikqkygnkiswadlmlltgnvalenmgfktlgfgggradtwqsdeavywgaettfvpqgndvrynnsvdinaradklekplaathmgliyvnpegpngtpdpaasakdireafgrmgmndtetvaliagghafgkthgavkgsnigpapeaadlgmqglgwhnsvgdgngpnqmtsgleviwtktptkwsngyleslinnnwtlvespagahqweavngtvdypdpfdktkfrkatmltsdlalindpeylkisqrwlehpeeladafakawfkllhrdlgpttrylgpevp"); // d3ut2a1
	s4 = sequence::from_string("lvhvasvekgrsyedfqkvynaialklreddeydnyigygpvlvrlawhisgtwdkhdntggsyggtyrfkkefndpsnaglqngfkflepihkefpwissgdlfslggvtavqemqgpkipwrcgrvdtpedttpdngrlpdadkdagyvrtffqrlnmndrevvalmgahalgkthlknsgyegpggaannvftnefylnllnedwklekndanneqwdsksgymmlptdysliqdpkylsivkeyandqdkffkdfskafekllengitfpkdapspfifktleeqgl"); // d2euta_

	benchmark_ungapped(s1, s2);
#ifdef __SSSE3__
	benchmark_ungapped_sse(s1, s2);
#endif
	benchmark_transpose();
#ifdef __SSE2__
	banded_swipe(s1, s2);
	swipe_cell_update();
	diag_scores(s1, s2);
#endif
#ifdef __SSE4_1__
	swipe(s3, s4);
#endif
}

}}