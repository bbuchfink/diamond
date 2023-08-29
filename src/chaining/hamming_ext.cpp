#include "chaining.h"
#include "../basic/config.h"
#include "../stats/score_matrix.h"
#include "../stats/stats.h"
#include "../util/util.h"

using std::stable_sort;
using std::vector;
using std::max;

namespace Chaining {

static ApproxHsp find_aln(vector<DiagonalSegment>::iterator begin, vector<DiagonalSegment>::iterator end, Loc qlen, Loc tlen) {
	for (auto it = begin; it < end; ++it) {
		const double ev = score_matrix.evalue(it->score, qlen, tlen);
		if ((it->id_percent() >= config.approx_min_id || Stats::approx_id(it->score, it->len, it->len) >= config.approx_min_id)
			&& max(it->cov_percent(qlen), it->cov_percent(tlen)) >= config.member_cover
			&& ev <= config.max_evalue)
			return ApproxHsp(0, 0, it->score, 0, it->query_range(), it->subject_range(), *it, ev);
	}
	return ApproxHsp(0);
}

static ApproxHsp filter(vector<DiagonalSegment>::iterator begin, vector<DiagonalSegment>::iterator end, Loc qlen, Loc tlen) {
	const double TOLERANCE_FACTOR = 1.1;
	stable_sort(begin, end, DiagonalSegment::cmp_score);
	Loc ident = 0, len = 0;
	const Loc qtol = safe_cast<Loc>(qlen * TOLERANCE_FACTOR), ttol = safe_cast<Loc>(tlen * TOLERANCE_FACTOR);
	for (auto it = begin; it < end; ++it) {
		if (len + it->len > qtol || len + it->len > ttol)
			continue;
		ident += it->ident;
		len += it->len;
	}
	if (max((double)len / qlen, (double)len / tlen) * 100.0 < config.diag_filter_cov
		|| (double)ident / len * 100.0 < config.diag_filter_id)
		return ApproxHsp(0, -1);
	return ApproxHsp(0);
}

ApproxHsp hamming_ext(vector<DiagonalSegment>::iterator begin, vector<DiagonalSegment>::iterator end, Loc qlen, Loc tlen) {
	if (config.hamming_ext) {
		ApproxHsp h = find_aln(begin, end, qlen, tlen);
		if (h.score > 0)
			return h;
	}
	if (config.diag_filter_cov > 0.0 || config.diag_filter_id > 0.0) {
		return filter(begin, end, qlen, tlen);
	}
	return ApproxHsp(0);
}

}