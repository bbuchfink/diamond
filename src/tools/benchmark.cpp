#include <chrono>
#include "../basic/sequence.h"
#include "../basic/score_matrix.h"
#include "../dp/score_vector.h"
#include "../util/simd/transpose.h"

using std::vector;
using std::chrono::high_resolution_clock;
using std::chrono::nanoseconds;
using std::chrono::duration_cast;

__m128i benchmark_global_128;

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
	static const size_t n = 100000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	const Letter *q = s1.data(), *s = s2.data();
	int score = 0;

	for (size_t i = 0; i < n; ++i) {

		score += xdrop_window2(q, s);

	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	std::chrono::nanoseconds time_span = duration_cast<std::chrono::nanoseconds>(t2 - t1);

	cout << score << endl;
	cout << "ns/Cell = " << (double)time_span.count() / (n*40) << endl;
}

void benchmark_ungapped_sse(const sequence &s1, const sequence &s2)
{
	static const size_t n = 100000000llu;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	const Letter *q = s1.data(), *s = s2.data();
	int score = 0;
	score_vector<uint8_t> sv, sum;

	for (size_t i = 0; i < n; ++i) {
		__m128i s = _mm_loadu_si128((const __m128i*)s1.data());
		sv.set_ssse3(i & 7, s);
		sum = sum + sv;
		sum.max(sv);
	}

	_mm_storeu_si128(&benchmark_global_128, sum.data_);
	cout << "ns/Cell = " << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 16) << endl;
}

void benchmark_transpose() {
	static const size_t n = 100000000llu;
	static char in[256], out[256];

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (size_t i = 0; i < n; ++i) {
		transpose(in, out, 0);
	}
	cout << "ns/Cell = " << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * 64) << endl;
}

void benchmark() {
	vector<Letter> s1, s2;
	
	s1 = sequence::from_string("QADATVATFFNGIDMPNQTNKTAAFLCAALGGPNAWTGRNLKEVHANMGVSNAQFTTVIGHLRSALTGAGVAAALVEQTVAVAETVRGDVVTV");
	s2 = sequence::from_string("QNDSSIIDFIKINDLAEQIEKISKKYIVSIVLGGGNIWRGSIAKELDMDRNLADNMGMMATIINGLALENALNHLNVNTIVLSAIKCDKLVHESSANNIKKAIEKEQVMIFVAGTGFPYFTTDSCAAIRAAETESSIILMGKNGVDGVYDSDPKINPNAQFYEHITFNMALTQNLKVMDATALALCQENNINLLVFNIDKPNAIVDVLEKKNKYTIVSK");
	
	benchmark_ungapped(s1, s2);
	benchmark_ungapped_sse(s1, s2);
	benchmark_transpose();
}