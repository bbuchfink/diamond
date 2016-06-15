#include "../basic/sequence.h"
#include "../basic/match.h"
#include "../basic/score_matrix.h"

Diagonal_segment score_diagonal(const Letter *query, const Letter *subject)
{
	int i = 0, j = 0, max_score = 0, score = 0, begin = 0, end = 0;
	while (query[i] != '\xff' && subject[i] != '\xff') {
		score += score_matrix(query[i], mask_critical(subject[i]));
		if (score <= 0) {
			score = 0;
			j = i + 1;
		}
		if (score > max_score) {
			max_score = score;
			begin = j;
			end = i + 1;
		}
		++i;
	}
	return Diagonal_segment(begin, begin, end - begin, max_score);
}

int score_range(sequence query, const Letter *subject, int i, int j, int j_end)
{
	int score = 0;
	while (j < j_end) {
		score += score_matrix(query[i], mask_critical(subject[j]));
		++i;
		++j;
	}
	return score;
}

void get_ungapped(sequence query, const Letter *subject, vector<Diagonal_segment> &out)
{
	const int min_score = 15;
	for (unsigned i = 0; i < query.length(); ++i) {
		Diagonal_segment d = score_diagonal(&query[i], subject);
		if (d.score >= min_score) {
			d.query_pos += i;
			out.push_back(d);
		}
	}
	unsigned i = 1;
	while (subject[i] != '\xff') {
		Diagonal_segment d = score_diagonal(&query[0], &subject[i]);
		if (d.score >= min_score) {
			d.subject_pos += i;
			out.push_back(d);
		}
		++i;
	}
}

int get_vgap_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, const Letter *subject)
{
	return 0;
}

int get_hgap_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, const Letter *subject)
{
	const int d = d2.diag() - d1.diag(), j2_end = d2.subject_pos + d2.len;
	int j1 = d1.subject_pos, j2 = j1 + d + 1, i1 = d1.query_pos, i2 = i1 + 1;
	int score = score_matrix(query[i1], mask_critical(subject[j1]))
		+ score_range(query, subject, i2, j2, j2_end)
		- config.gap_open - d*config.gap_extend;
	int max_score = 0;
	while (true) {
		max_score = std::max(score, max_score);
		score -= score_matrix(query[i2], mask_critical(subject[j2]));
		++i1; ++i2; ++j1; ++j2;
		if (j2 >= j2_end)
			break;
		score += score_matrix(query[i1], mask_critical(subject[j1]));
	}
	return max_score;
}

int get_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, const Letter *subject)
{
	if (d1.diag() > d2.diag())
		return get_vgap_link(d1, d2, query, subject);
	else
		return get_hgap_link(d1, d2, query, subject);
}

void get_links(const vector<Diagonal_segment> &diag, sequence query, const Letter *subject)
{
	for (vector<Diagonal_segment>::const_iterator i = diag.begin(); i != diag.end();++i)
		for (vector<Diagonal_segment>::const_iterator j = diag.begin(); j != diag.end(); ++j) {
			int link_score = get_link(*i, *j, query, subject);
#ifdef LOG_GA
			if (link_score > i->score && link_score > j->score)
				cout << i - diag.begin() << ' ' << j - diag.begin() << ' ' << link_score << endl;
#endif
		}
}

void greedy_align(sequence query, local_match &hsp)
{
	vector<Diagonal_segment> diag;
	get_ungapped(query, hsp.subject_, diag);
#ifdef LOG_GA
	for (vector<Diagonal_segment>::const_iterator i = diag.begin(); i != diag.end(); ++i)
		cout << i-diag.begin() << ' ' << i->query_pos << ' ' << i->subject_pos << ' ' << i->score << endl;
	cout << endl;
#endif
	get_links(diag, query, hsp.subject_);
}