#include <chrono>
#include "../dp/score_vector_int8.h"
#include "../dp/score_vector_int16.h"
#include "../dp/swipe/cell_update.h"

using namespace DISPATCH_ARCH;
using std::vector;
using std::chrono::high_resolution_clock;
using std::chrono::nanoseconds;
using std::chrono::duration_cast;
using std::cout;
using std::endl;
using std::array;

namespace Benchmark { namespace DISPATCH_ARCH {

#ifdef __SSE4_1__

static const size_t N = 256;
using Sv = ScoreVector<int8_t, SCHAR_MIN>;
using Cell = ForwardCell<Sv>;
//using Cell = Sv;
static char query[N];
static const size_t C = ScoreTraits<Sv>::CHANNELS;

static void update_row(Cell* diagonal_cell, Cell* horizontal_gap, Sv* profile) {
	//auto id_mask = DummyIdMask<Sv>(0, Sv());
	auto id_mask = VectorIdMask<Sv>(0, Sv());
	//auto row_counter = DummyRowCounter(0);
	auto row_counter = VectorRowCounter<Sv>(0);
	Sv scores, gap_extension, gap_open, best;
	Cell vertical_gap;

	for (size_t i = 0; i < N; ++i) {
		diagonal_cell[i] = swipe_cell_update(diagonal_cell[i], profile[(int)query[i]], nullptr, gap_extension, gap_open, horizontal_gap[i], vertical_gap, best, nullptr, row_counter, id_mask);
	}
}
	
void swipe_cell_update() {

	static const size_t n = 1000000llu;
	Cell diagonal_cell[N], horizontal_gap[N];
	Sv profile[32];
	for (size_t i = 0; i < N; ++i) {
		query[i] = rand() % 32;
	}
	for (size_t i = 0; i < 32; ++i) {
		int8_t v[C];
		for (size_t j = 0; j < C; ++j)
			v[j] = rand() % 20 - 10;
		profile[i] = load_sv<Sv>(v);
	}
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	for (size_t i = 0; i < n; ++i) {
		update_row(diagonal_cell, horizontal_gap, profile);
	}
	volatile auto x = diagonal_cell[0].data_;
	volatile auto y = diagonal_cell[0].ident;
	volatile auto z = diagonal_cell[0].len;

	cout << "SWIPE cell update (int8_t):\t" << (double)duration_cast<std::chrono::nanoseconds>(high_resolution_clock::now() - t1).count() / (n * N * 32) * 1000 << " ps/Cell" << endl;
}

#endif

}}