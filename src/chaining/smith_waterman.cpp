#include "chaining.h"
#include "../output/output_format.h"
#include "../dp/dp.h"

void print_diag(int i0, int j0, int l, int score, const Diag_graph &diags, const sequence &query, const sequence &subject)
{
	Diagonal_segment ds(i0, j0, l, 0);
	unsigned n = 0;
	int path_max, path_min;
	for (vector<Diagonal_node>::const_iterator d = diags.nodes.begin(); d != diags.nodes.end(); ++d) {
		if (d->intersect(ds).len > 0) {
			if (d->score == 0)
				continue;
			const int diff = score_range(query, subject, d->query_end(), d->subject_end(), j0 + l);
			if (n > 0)
				cout << "(";
			cout << "Diag n=" << d - diags.nodes.begin() << " i=" << i0 << " j=" << j0 << " len=" << l
				<< " prefix_score=" << score + score_range(query, subject, i0 + l, j0 + l, d->subject_end()) - std::min(diff, 0)
				<< " prefix_score2=" << diags.prefix_score((unsigned)(d - diags.nodes.begin()), j0 + l, path_max, path_min);
			if (n > 0)
				cout << ")";
			cout << endl;
			++n;
		}
	}
	if (n == 0)
		cout << "Diag n=x i=" << i0 << " j=" << j0 << " len=" << l << " prefix_score=" << score << endl;
}

void smith_waterman(sequence q, sequence s, const Diag_graph &diags)
{
	Hsp hsp;
	smith_waterman(q, s, hsp);
	Hsp::Iterator i = hsp.begin();
	int i0 = -1, j0 = -1, l = 0, score = 0;
	for (; i.good(); ++i) {
		switch (i.op()) {
		case op_match:
		case op_substitution:
			if (i0 < 0) {
				i0 = i.query_pos.translated;
				j0 = i.subject_pos;
				l = 0;
			}
			score += score_matrix(q[i.query_pos.translated], s[i.subject_pos]);
			++l;
			break;
		case op_deletion:
		case op_insertion:
			if (i0 >= 0) {
				print_diag(i0, j0, l, score, diags, q, s);
				score -= score_matrix.gap_open() + score_matrix.gap_extend();
				i0 = -1;
				j0 = -1;
			}
			else
				score -= score_matrix.gap_extend();
			break;
		case op_frameshift_forward:
		case op_frameshift_reverse:
			;
		}
	}
	print_diag(i0, j0, l, score, diags, q, s);
	print_hsp(hsp, TranslatedSequence(q));
}