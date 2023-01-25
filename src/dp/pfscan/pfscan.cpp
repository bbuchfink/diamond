#include "pfscan.h"
#include "smith_waterman.h"
#include "../../util/sequence/sequence.h"
#include "../ungapped.h"

using std::min;

namespace DP { namespace PrefixScan { namespace DISPATCH_ARCH {

Hsp align16(const Config& cfg) {
	using SC = StaticConfig<::DISPATCH_ARCH::ScoreVector<int16_t, 0>, Local>;
	Hsp hsp = banded_smith_waterman<SC>(cfg);
	//hsp.evalue = score_matrix.evalue(hsp.score, cfg.query.length(), cfg.target.length());
	return hsp;
}

Hsp align8(const Config& cfg) {
#ifdef __SSE4_1__
	using SC = StaticConfig<::DISPATCH_ARCH::ScoreVector<int8_t, 0>, Local>;
	Hsp hsp = banded_smith_waterman<SC>(cfg);
	//hsp.evalue = score_matrix.evalue(hsp.score, cfg.query.length(), cfg.target.length());
	return hsp;
#else
	return Hsp();
#endif
}

static Hsp align_dispatch_score_16(const Config& cfg) {
	using SC = StaticConfig<::DISPATCH_ARCH::ScoreVector<int16_t, 0>, Anchored>;
	Config cfg16(cfg);
	cfg16.adjust_band(16);
	cfg16.score_bias = score_matrix.gap_extend() * (cfg16.d_end - cfg16.d_begin) * 2;
	const Hsp h = banded_smith_waterman<SC>(cfg16);
	cfg.stats.inc(Statistics::EXT16);
	if (h.score < SCHAR_MAX)
		cfg.stats.inc(Statistics::EXT_WASTED_16);
	if (h.score == (Score)SHRT_MAX - cfg16.score_bias) {
		task_timer timer;
		using SC = StaticConfig<int, Anchored>;
		const Hsp h = banded_smith_waterman<SC>(cfg);
		cfg.stats.inc(Statistics::TIME_EXT_32, timer.microseconds());
		cfg.stats.inc(Statistics::EXT32);
		return h;
	}
	return h;
}

static Hsp align_dispatch_score(const Config& cfg) {
	if (cfg.query.length() >= SHRT_MAX || cfg.target.length() >= SHRT_MAX) {
		task_timer timer;
		using SC = StaticConfig<int, Anchored>;
		const Hsp h = banded_smith_waterman<SC>(cfg);
		cfg.stats.inc(Statistics::TIME_EXT_32, timer.microseconds());
		cfg.stats.inc(Statistics::EXT32);
		return h;
	}
	else {
#ifdef __SSE4_1__
		if (config.no_8bit_extension || cfg.band() / 32 > SCHAR_MAX || cfg.band() <= 16 || cfg.score_hint >= 95 || cfg.band() > 128
			|| score_matrix.gap_open() + score_matrix.gap_extend() * max(-cfg.d_begin, cfg.d_end) > SCHAR_MAX) {
#endif
			return align_dispatch_score_16(cfg);
#ifdef __SSE4_1__
		}
		else {
			using SC = StaticConfig<::DISPATCH_ARCH::ScoreVector<int8_t, 0>, Anchored>;
			Config cfg8(cfg);
			cfg8.adjust_band(32);
			const Hsp h = banded_smith_waterman<SC>(cfg8);
			cfg.stats.inc(Statistics::EXT8);
			if (h.score == (Score)SCHAR_MAX) {
				cfg.stats.inc(Statistics::EXT_OVERFLOW_8);
				return align_dispatch_score_16(cfg);
			}
			return h;
		}
#endif
	}
}

static double length_fraction(Loc pos, Loc len, Loc anchor_len) {
	return double(len - pos) / (len - anchor_len);
}

static Hsp align_right(Loc i, Loc j, Loc d_begin, Loc d_end, Score prefix_score, const Config& cfg) {
	const Sequence query(cfg.query.subseq(i, cfg.query.length())),
		target(cfg.target.subseq(j, cfg.target.length()));
	vector<const int16_t*> profile;
	vector<const int8_t*> profile8;
	profile.reserve(AMINO_ACID_COUNT);
	profile8.reserve(AMINO_ACID_COUNT);
	for (size_t a = 0; a < AMINO_ACID_COUNT; ++a) {
		profile.push_back(cfg.query_profile[a] + i);
		profile8.push_back(cfg.query_profile8[a] + i);
	}
	const int band = std::max(32, Loc((d_end - d_begin) * 0.15));
	d_begin -= band;
	d_end += band - 1;
	const Loc d0 = Geo::clip_diag(Geo::diag_sub_matrix(d_begin, i, j), query.length(), target.length()),
		d1 = Geo::clip_diag(Geo::diag_sub_matrix(d_end, i, j), query.length(), target.length());
	Config cfg_r{
		query,
		target,
		cfg.query_seqid,
		cfg.target_seqid,
		d0,
		d1,
		cfg.hint_target_range.length() > 0 ? Interval(0, cfg.hint_target_range.end_ - j) : Interval(),
		profile.data(),
		nullptr,
		profile8.data(),
		nullptr,
		cfg.stats,
		0,
		std::max(score_matrix.gap_extend() * (d1 - d0), score_matrix.rawscore(config.gapped_xdrop)),
		//score_matrix.rawscore(config.gapped_xdrop)
		cfg.score_hint - prefix_score
	};
	Hsp h = align_dispatch_score(cfg_r);
	h.query_range.end_ += i;
	h.subject_range.end_ += j;
	return h;
}

static Hsp align_left(Loc i, Loc j, Loc d_begin, Loc d_end, Score suffix_score, const Config& cfg) {
	const vector<Letter> ql(cfg.query.reverse()),
		tl(cfg.target.reverse());
	const Loc qlen = cfg.query.length(), tlen = cfg.target.length();
	Config cfg_l{
		Sequence(ql),
		Sequence(tl),
		cfg.query_seqid,
		cfg.target_seqid,
		cfg.d_begin,
		cfg.d_end,		
		cfg.hint_target_range.length() > 0 ? Interval(0, tlen - cfg.hint_target_range.begin_) : Interval(),
		cfg.query_profile_rev,
		nullptr,
		cfg.query_profile_rev8,
		nullptr,
		cfg.stats,
		cfg.score_bias,
		cfg.xdrop,
		cfg.score_hint
	};
	Hsp h = align_right(qlen - 1 - i, tlen - 1 - j, Geo::rev_diag(d_end - 1, qlen, tlen), Geo::rev_diag(d_begin, qlen, tlen) + 1, suffix_score, cfg_l);
	h.query_range.begin_ = qlen - 1 - (h.query_range.end_ - 1);
	h.subject_range.begin_ = tlen - 1 - (h.subject_range.end_ - 1);
	return h;
}

Hsp align_anchored(const Anchor& anchor, const Config& cfg) {
	assert(anchor.diag() >= cfg.d_begin && anchor.diag() < cfg.d_end);
	Hsp h(false, anchor.score);
	h.query_range = anchor.query_range();
	h.subject_range = anchor.subject_range();
	if (anchor.query_end() < cfg.query.length() && anchor.subject_end() < cfg.target.length()) {
		Hsp r = align_right(anchor.query_end(), anchor.subject_end(), anchor.d_min_right, anchor.d_max_right + 1, anchor.prefix_score, cfg);
		h.score += r.score;
		h.query_range.end_ = r.query_range.end_;
		h.subject_range.end_ = r.subject_range.end_;
	}
	if (anchor.query_begin() > 0 && anchor.subject_begin() > 0) {
		const Score suffix_score = cfg.score_hint - anchor.prefix_score + anchor.score;
		Hsp l = align_left(anchor.query_begin() - 1, anchor.subject_begin() - 1, anchor.d_min_left, anchor.d_max_left + 1, suffix_score, cfg);
		h.score += l.score;
		h.query_range.begin_ = l.query_range.begin_;
		h.subject_range.begin_ = l.subject_range.begin_;
	}
	h.query_source_range = h.query_range;
	h.bit_score = score_matrix.bitscore(h.score);
	h.evalue = score_matrix.evalue(h.score, cfg.query.length(), cfg.target.length());
	if (h.evalue > config.max_evalue) {
		h.evalue = DBL_MAX;
		h.score = 0;
	}
	return h;
}

}}}