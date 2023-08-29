#include "pfscan.h"
#include "../../basic/config.h"
#include "../../util/geo/geo.h"
#include "../dp/score_vector_int8.h"
#include "../dp/score_vector_int16.h"
#include "../../util/memory/alignment.h"
#include "simd.h"

using std::max;
using std::min;
using std::vector;
using std::fill;
using std::numeric_limits;
using std::tie;

namespace DP { namespace PrefixScan { namespace DISPATCH_ARCH {

static const double XDROP_MAX_MASKED_RATIO = 0.5;

struct Local {};
struct Anchored {};

template<typename ScoreVector, typename L>
struct StaticConfig {
	using Sv = ScoreVector;
	using Logic = L;
};

template<typename Sv>
static inline void saturate(Sv& sv, const Sv& zero, Anchored) {
}

template<typename Sv>
static inline void saturate(Sv& sv, const Sv& zero, Local) {
	sv = max(sv, zero);
}

template<typename Score>
static inline void init_scores(Score* score, Loc d, Loc band, Score bias, Local) {
	fill(score, score + band, (Score)bias);
}

template<typename Score>
static inline void init_scores(Score* score, Loc d, Loc band, Score bias, Anchored) {
	const Score e = score_matrix.gap_extend();
	Score s = bias - score_matrix.gap_open();
	for (Loc i = d - 1; i >= 0; --i) {
		s -= e;
		score[i] = s;
	}
	score[d] = bias;
	s = bias - score_matrix.gap_open();
	for (Loc i = d + 1; i < band; ++i) {
		s -= e;
		score[i] = s;
	}
}

template<typename Score>
static inline const int16_t* const* get_profile(const Config& cfg, Score) {
	return cfg.query_profile;
}

static inline const int8_t* const* get_profile(const Config& cfg, int8_t) {
	return cfg.query_profile8;
}

template<typename Sv, typename Score>
static inline Sv load_profile(const Score* ptr) {
	return Sv(ptr);
}

template<>
inline int32_t load_profile(const int16_t* ptr) {
	return *ptr;
}

template<Loc channels>
static inline Loc clip_i0(Loc i0) {
	return max(i0 + (-i0 / channels) * channels, i0);
}

Statistics::value cell_stat(int8_t) {
	return Statistics::DP_CELLS_8;
}

Statistics::value cell_stat(int16_t) {
	return Statistics::DP_CELLS_16;
}

Statistics::value cell_stat(int32_t) {
	return Statistics::DP_CELLS_32;
}

template<typename StaticConfig>
Hsp FLATTEN banded_smith_waterman(const Config& cfg) {
	using Sv = typename StaticConfig::Sv;
	using Score = typename ::DISPATCH_ARCH::ScoreTraits<Sv>::Score;
	using Logic = typename StaticConfig::Logic;

	constexpr Loc SCORE_MIN = (Loc)numeric_limits<Score>::min(), SCORE_MAX = (Loc)numeric_limits<Score>::max();
	const Loc channels = ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS;

	TaskTimer timer;

	const Loc band = cfg.d_end - cfg.d_begin,
		qlen = cfg.query.length(),
		tlen = cfg.target.length(),
		j0 = max(Geo::j(0, cfg.d_end - 1), 0),
		j1 = min(Geo::j(qlen, cfg.d_begin), tlen);
	Loc i0 = Geo::i(j0, cfg.d_begin),
		i1 = Geo::i(j0, cfg.d_end);

	assert(band % channels == 0);
	auto query_profile = get_profile(cfg, Score());
	vector<Score, Util::Memory::AlignmentAllocator<Score, 32>> scores(band);
	vector<Score, Util::Memory::AlignmentAllocator<Score, 32>> hgap(band + 1, numeric_limits<Score>::min() + score_matrix.gap_extend());
	init_scores(scores.data(), -i0, band, (Score)cfg.score_bias, Logic());
	const Sv gap_extend = Sv(-score_matrix.gap_extend()),
		gap_open = Sv(-score_matrix.gap_open()),
		zero(cfg.score_bias);
	Sv pf_const2, pf_const1, max_score((Score)cfg.score_bias), max_j(-1), max_i(-1);
	tie(pf_const1, pf_const2) = prefix_scan_consts(gap_extend);

	Loc j = j0;
	size_t n = 0;
	for (; j < j1; ++j) {
		const Loc i0_ = clip_i0<channels>(i0);
		const auto sv_offset = (i0_ - i0);
		Score* score_it = scores.data() + sv_offset;
		Score* hgap_it = hgap.data() + sv_offset;
		auto profile = query_profile[cfg.target[j]] + i0_;
		Sv vgap(numeric_limits<Score>::min() + score_matrix.gap_extend());
		Sv col_max(SHRT_MIN), col_max_i(-1), counter(0);
		for (Loc i = i0_; i < i1; i += channels) {
			Sv score = load_sv_aligned<Sv>(score_it) + load_profile<Sv>(profile);
			Sv hgap = load_sv<Sv>(hgap_it + 1);
			hgap += gap_extend;
			score = max(score, hgap);

			Sv v = score + gap_open;
			v = prefix_scan(v, gap_extend, pf_const2);
			v = max(v, vgap + pf_const1);
			score = max(score, v);
			saturate(score, zero, Logic());

			store_aligned(score, score_it);
			
			vgap = Sv(extract<channels - 1>(v));

			hgap = max(hgap, score + gap_open);
			store_aligned(hgap, hgap_it);

			const Sv gt_mask = score > col_max;
			col_max_i = blend(col_max_i, counter, gt_mask);
			col_max = max(col_max, score);

			score_it += channels;
			hgap_it += channels;
			profile += channels;
			counter += Sv(1);
			++n;
		}
		++i0;
		i1 = min(i1 + 1, qlen);
		const Sv gt_mask = col_max > max_score;
		max_j = blend(max_j, Sv(min(j + SCORE_MIN, SCORE_MAX)), gt_mask);
		max_i = blend(max_i, col_max_i, gt_mask);
		max_score = max(max_score, col_max);

		if ((j & 31) == 31) {
			const auto best = max_entry(max_score);
			if (best.first == numeric_limits<Score>::max())
				break;
			if (j < cfg.hint_target_range.end_)
				continue;
			const auto col_best = max_entry(col_max);
			if ((int32_t)best.first - (int32_t)col_best.first >= cfg.xdrop
				&& cfg.target.subseq_clipped(j - cfg.xdrop, j + 1).masked_letter_ratio() < XDROP_MAX_MASKED_RATIO
				&& cfg.query.subseq_clipped(i1 - cfg.xdrop, i1 + 1).masked_letter_ratio() < XDROP_MAX_MASKED_RATIO)
				break;
		}
	}
	Hsp out;
	const auto best = max_entry(max_score);
	out.score = ::DISPATCH_ARCH::ScoreTraits<Sv>::int_score(best.first) - cfg.score_bias;
	if (out.score > 0) {
		const Loc max_j_ = (Loc)extract(max_j, best.second);
		if (max_j_ == SCORE_MAX)
			out.score = SCORE_MAX;
		out.subject_range.end_ = max_j_ - SCORE_MIN + 1;
		out.query_range.end_ = clip_i0<channels>(i0 - j + out.subject_range.end_ - 1) + channels * extract(max_i, best.second) + best.second + 1;
		assert(out.query_range.end_ > 0 && out.subject_range.end_ > 0);
	}
	cfg.stats.inc(cell_stat(Score()), n * channels);
	cfg.stats.inc(Statistics::TIME_SW, timer.microseconds());
	return out;
}

}}}