#pragma once
#include "../basic/sequence.h"
#include "../stats/score_matrix.h"
#include "../util/util.h"

inline int nw_semiglobal(Sequence query, Sequence target) {
	Matrix<int> m(query.length() + 1, target.length() + 1);
	std::vector<int> hgap(query.length() + 1);
	for (Loc i = 1; i <= query.length(); ++i) {
		hgap[i] = -score_matrix.gap_open() - i * score_matrix.gap_extend();
		m[i][0] = hgap[i] - score_matrix.gap_open() - score_matrix.gap_extend();
	}
	int score = INT_MIN;
	for (int j = 1; j <= target.length(); ++j) {
		int vgap = -score_matrix.gap_open() - score_matrix.gap_extend();
		for (int i = 1; i <= query.length(); ++i) {
			int s = m[i - 1][j - 1] + score_matrix(query[i - 1], target[j - 1]);
			s = std::max(s, vgap);
			s = std::max(s, hgap[i]);
			m[i][j] = s;
			vgap -= score_matrix.gap_extend();
			hgap[i] -= score_matrix.gap_extend();
			int open = s - score_matrix.gap_open() - score_matrix.gap_extend();
			vgap = std::max(vgap, open);
			hgap[i] = std::max(hgap[i], open);
		}
		score = std::max(score, m[query.length()][j]);
	}
	return score;
}