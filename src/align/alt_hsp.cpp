#include "target.h"
#include "../dp/dp.h"
#include "../util/sequence/sequence.h"

using std::list;
using std::array;
using std::vector;

namespace Extension {

struct ActiveTarget {
	ActiveTarget(vector<Match>::iterator match, SequenceSet& dst) :
		match(match),
		active(0)
	{
		masked_seq.fill(nullptr);
		uint32_t active = 0;
		const auto l = match->seq.length();
		for (const Hsp& h : match->hsp) {
			if (!(active & (1 << h.frame))) {
				dst.reserve(l);
				active |= 1 << h.frame;
			}
		}
	}
	ActiveTarget(const ActiveTarget& t):
		match(t.match),
		masked_seq(t.masked_seq),
		active(0)
	{
		for (int32_t i = 0; i < align_mode.query_contexts; ++i)
			if (!(t.active & (1 << i)))
				masked_seq[i] = nullptr;
	}
	void copy_seq(SequenceSet& dst, int64_t& i) {
		for (const Hsp& h : match->hsp) {
			if (!masked_seq[h.frame]) {
				dst.assign(i, match->seq.data(), match->seq.end());
				masked_seq[h.frame] = dst.ptr(i++);
			}
			Letter* seq = masked_seq[h.frame];
			std::fill(seq + h.subject_range.begin_, seq + h.subject_range.end_, SUPER_HARD_MASK);
		}
	}
	Sequence masked(int32_t context) const {
		return Sequence(masked_seq[context], match->seq.length());
	}
	int32_t check_fully_masked() {
		int32_t n = 0;
		for (int32_t i = 0; i < align_mode.query_contexts; ++i)
			if (active & (1 << i)) {
				if (Util::Seq::is_fully_masked(masked(i)))
					active &= ~(1 << i);
				else
					++n;
			}
		return n;
	}
	const vector<Match>::iterator match;
	array<Letter*, MAX_CONTEXT> masked_seq;
	uint32_t active;
};

using TargetVec = vector<ActiveTarget>;

static TargetVec recompute_alt_hsps(const Sequence* query_seq, const int query_source_len, const Bias_correction* query_cb, TargetVec& targets, const HspValues v, Statistics& stats) {
	array<DP::Targets, MAX_CONTEXT> dp_targets;
	const Loc qlen = query_seq[0].length();
	for (auto it = targets.begin(); it != targets.end(); ++it) {
		const int64_t dp_size = (int64_t)qlen * (int64_t)it->match->seq.length();
		const ::Stats::TargetMatrix* matrix = it->match->matrix.blank() ? nullptr : &it->match->matrix;
		const int bin = DP::BandedSwipe::bin(v, qlen, 0, 0, dp_size, matrix ? matrix->score_width() : 0, 0);
		for (int32_t context = 0; context < align_mode.query_contexts; ++context) {
			if (it->masked_seq[context]) {
				const Sequence target = it->masked(context);
				dp_targets[context][bin].emplace_back(target, target.length(), BlockId(it - targets.begin()), matrix);
			}
		}
	}

	for (int32_t context = 0; context < align_mode.query_contexts; ++context) {
		const int8_t* cbs = ::Stats::CBS::hauser(config.comp_based_stats) ? query_cb[context].int8.data() : nullptr;
		DP::Params params{ query_seq[context], "", Frame(context), query_source_len, cbs, DP::Flags::FULL_MATRIX, v, stats, nullptr };
		list<Hsp> hsp = DP::BandedSwipe::swipe(dp_targets[context], params);
		while (!hsp.empty()) {
			ActiveTarget& t = targets[hsp.front().swipe_target];
			list<Hsp>& l = t.match->hsp;
			l.splice(l.end(), hsp, hsp.begin());
			std::fill(t.masked_seq[context] + l.back().subject_range.begin_, t.masked_seq[context] + l.back().subject_range.end_, SUPER_HARD_MASK);
			t.active |= 1 << context;
		}
	}

	TargetVec out;
	for (ActiveTarget& t : targets) {
		if (t.active) {
			t.match->inner_culling();
			if (t.check_fully_masked() > 0 && (t.match->hsp.size() < config.max_hsps || config.max_hsps == 0))
				out.emplace_back(t);
		}
	}
	return out;
}

void recompute_alt_hsps(vector<Match>::iterator begin, vector<Match>::iterator end, const Sequence* query, const int query_source_len, const Bias_correction* query_cb, const HspValues v, Statistics& stats) {
	if (config.max_hsps == 1)
		return;
	TargetVec targets;
	targets.reserve(end - begin);
	SequenceSet target_seqs;
	for (auto i = begin; i != end; ++i)
		targets.emplace_back(i, target_seqs);
	target_seqs.finish_reserve();
	int64_t i = 0;
	for (ActiveTarget& t : targets)
		t.copy_seq(target_seqs, i);

	while(!targets.empty())
		targets = recompute_alt_hsps(query, query_source_len, query_cb, targets, v, stats);
}

}