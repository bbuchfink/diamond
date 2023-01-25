#include "score_profile.h"
#include "score_vector.h"
#include "score_vector_int8.h"
#include "score_vector_int16.h"
#include "../util/util.h"

using std::array;

namespace DP { namespace DISPATCH_ARCH {

template<typename Score>
LongScoreProfile<Score> make_profile(Sequence seq, const int8_t* cbs, int64_t padding)
{
	LongScoreProfile<Score> p;
	p.padding = std::max(padding, (int64_t)LongScoreProfile<Score>::DEFAULT_PADDING);
	for (unsigned l = 0; l < AMINO_ACID_COUNT; ++l) {
		const int8_t* scores = &score_matrix.matrix8()[l << 5];
		p.data[l].reserve(round_up(seq.length(), 32) + 2 * p.padding);
		p.data[l].insert(p.data[l].end(), p.padding, -1);
#if ARCH_ID == 2
		using Sv = ::DISPATCH_ARCH::ScoreVector<int8_t, 0>;
		constexpr auto CHANNELS = ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS;
		alignas(32) array<Score, CHANNELS> buf;
		for (Loc i = 0; i < seq.length(); i += CHANNELS) {
			__m256i s = _mm256_loadu_si256((const __m256i*)(seq.data() + i));
			Sv scores(l, s);
			if (cbs && l < TRUE_AA)
				scores += Sv(cbs + i);
			store_expanded(scores, buf.data());
			p.data[l].insert(p.data[l].end(), buf.begin(), buf.end());
		}
		p.data[l].erase(p.data[l].end() - round_up(seq.length(), (Loc)CHANNELS) + seq.length(), p.data[l].end());
#else
		for (Loc i = 0; i < seq.length(); ++i) {
			int8_t score = scores[(int)seq[i]];
			if (cbs)
				score += cbs[i];
			p.data[l].push_back(score);
		}
#endif
		p.data[l].insert(p.data[l].end(), p.padding, -1);
	}
	return p;
}

LongScoreProfile<int8_t> make_profile8(Sequence seq, const int8_t* cbs, int64_t padding) {
	return make_profile<int8_t>(seq, cbs, padding);
}

LongScoreProfile<int16_t> make_profile16(Sequence seq, const int8_t* cbs, int64_t padding) {
	return make_profile<int16_t>(seq, cbs, padding);
}

}}