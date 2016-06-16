#include <list>
#include "../basic/sequence.h"
#include "../basic/match.h"
#include "../basic/score_matrix.h"

#define LOG_GA

using std::list;

struct Link
{
	unsigned target;
	int subject_pos;
};

typedef vector<Link> Link_list;

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

int score_range(sequence query, sequence subject, int i, int j, int j_end)
{
	int score = 0;
	while (j < j_end) {
		score += score_matrix(query[i], mask_critical(subject[j]));
		++i;
		++j;
	}
	return score;
}

void get_ungapped(sequence query, sequence subject, vector<Diagonal_segment> &out)
{
	const int min_score = 15;
	for (unsigned i = 0; i < query.length(); ++i) {
		Diagonal_segment d = score_diagonal(&query[i], &subject[0]);
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

int get_vgap_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, sequence subject, Link &l)
{
	return 0;
}

int get_hgap_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query,  sequence subject, Link &l)
{
	const int d = d2.diag() - d1.diag(), j2_end = d2.subject_pos + d2.len;
	int j1 = d1.subject_pos, j2 = j1 + d + 1, i1 = d1.query_pos, i2 = i1 + 1;
	int score = score_matrix(query[i1], mask_critical(subject[j1]))
		+ score_range(query, subject, i2, j2, j2_end)
		- config.gap_open - d*config.gap_extend;
	int max_score = 0;
	while (true) {
		if (score > max_score) {
			max_score = score;
			l.subject_pos = j1;
		}		
		score -= score_matrix(query[i2], mask_critical(subject[j2]));
		++i1; ++i2; ++j1; ++j2;
		if (j2 >= j2_end)
			break;
		score += score_matrix(query[i1], mask_critical(subject[j1]));
	}
	return max_score;
}

int get_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, sequence subject, Link &l)
{
	if (d1.diag() > d2.diag())
		return get_vgap_link(d1, d2, query, subject, l);
	else
		return get_hgap_link(d1, d2, query, subject, l);
}

void get_links(const vector<Diagonal_segment> &diag, sequence query, sequence subject, vector<Link_list> &links, vector<bool> &is_root)
{
	for (vector<Diagonal_segment>::const_iterator i = diag.begin(); i != diag.end();++i)
		for (vector<Diagonal_segment>::const_iterator j = diag.begin(); j != diag.end(); ++j) {
			Link l;
			int link_score = get_link(*i, *j, query, subject, l);
			if (link_score > i->score && link_score > j->score) {
				l.target = (unsigned)(j - diag.begin());
				links[i - diag.begin()].push_back(l);
				is_root[l.target] = false;
#ifdef LOG_GA
				cout << i - diag.begin() << ' ' << j - diag.begin() << ' ' << link_score << endl;
#endif
			}
		}
}

std::ostream& indent(std::ostream &str, unsigned n)
{
	for (unsigned i = 0; i < n; ++i)
		str << ' ';
	return str;
}

void follow_path(unsigned level, unsigned node, vector<Link_list> &links, int score, int subject_pos, sequence query, sequence subject, const vector<Diagonal_segment> &diag)
{
#ifdef LOG_GA
	indent(cout, level) << "Node " << node << " Score=" << score << endl;
#endif
	const Link_list &l = links[node];
	const Diagonal_segment& d = diag[node];
	const int diff = subject_pos - d.subject_pos;
	if (l.size() == 0) {
		score += score_range(query, subject, d.query_pos + diff, subject_pos, d.subject_pos + d.len);
		indent(cout, level) << "Final score=" << score << endl << endl;
	}
	for (Link_list::const_iterator i = l.begin(); i != l.end(); ++i) {
		if (i->subject_pos < subject_pos)
			continue;
		const int shift = diag[i->target].diag() - d.diag();
		follow_path(level + 1,
			i->target,
			links,
			score + score_range(query, subject, d.query_pos + diff, subject_pos, i->subject_pos + 1),
			shift < 0 ? subject_pos + 1 : subject_pos + shift + 1,
			query,
			subject,
			diag);
	}
}

void greedy_align(sequence query, sequence subject)
{
	vector<Diagonal_segment> diag;
	vector<Link_list> links;
	get_ungapped(query, subject, diag);
#ifdef LOG_GA
	cout << "Diagonals:" << endl;
	for (vector<Diagonal_segment>::const_iterator i = diag.begin(); i != diag.end(); ++i)
		cout << i-diag.begin() << ' ' << i->query_pos << ' ' << i->subject_pos << ' ' << i->score << endl;
	cout << endl;
#endif
	links.resize(diag.size());
	vector<bool> is_root(diag.size(), true);
#ifdef LOG_GA
	cout << "Links:" << endl;
#endif
	get_links(diag, query, subject, links, is_root);
	for (unsigned i = 0; i < diag.size(); ++i)
		if (is_root[i])
			follow_path(0, i, links, 0, diag[i].subject_pos, query, subject, diag);
}