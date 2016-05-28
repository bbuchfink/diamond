#include "dp.h"
#include "../basic/score_matrix.h"

int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned seed_len, unsigned &delta, unsigned &len)
{
	int score(0), st(0);
	unsigned n(0);
	delta = 0;

	const Letter *q(query - 1), *s(subject - 1);
	const unsigned window_left = std::max(config.window, (unsigned)Const::seed_anchor) - Const::seed_anchor;
	while (score - st < config.xdrop
		&& delta < window_left
		&& *q != '\xff'
		&& *s != '\xff')
	{
		st += score_matrix(*q, mask_critical(*s));
		score = std::max(score, st);
		--q;
		--s;
		++delta;
	}

	q = query + seed_len;
	s = subject + seed_len;
	st = score;
	assert(seed_len >= Const::seed_anchor);
	const unsigned window_right = std::max(config.window, seed_len - Const::seed_anchor) - (seed_len - Const::seed_anchor);
	while (score - st < config.xdrop
		&& n < window_right
		&& *q != '\xff'
		&& *s != '\xff')
	{
		st += score_matrix(*q, mask_critical(*s));
		score = std::max(score, st);
		++q;
		++s;
		++n;
	}

	for (unsigned i = 0; i<seed_len; ++i)
		score += score_matrix(query[i], mask_critical(subject[i]));

	len = delta + n + seed_len;
	return score;
}

void xdrop_ungapped(const Letter *query, const Letter *subject, unsigned &delta, unsigned &len)
{
	int score(0), st(0), n=1;
	delta = 0;
	
	const Letter *q(query - 1), *s(subject - 1);
	while (score - st < config.xdrop
		&& *q != '\xff'
		&& *s != '\xff')
	{
		st += score_matrix(*q, mask_critical(*s));
		if (st > score) {
			score = st;
			delta = n;
		}
		--q;
		--s;
		++n;
	}

	q = query;
	s = subject;
	st = score;
	n = 1;
	len = 0;
	while (score - st < config.xdrop
		&& *q != '\xff'
		&& *s != '\xff')
	{
		st += score_matrix(*q, mask_critical(*s));
		if (st > score) {
			score = st;
			len = n;
		}
		++q;
		++s;
		++n;
	}
	len += delta;
}