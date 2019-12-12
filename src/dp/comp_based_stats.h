#ifndef COMP_BASED_STATS_H_
#define COMP_BASED_STATS_H_

#include <vector>
#include <stdint.h>
#include "../basic/sequence.h"
#include "../basic/diagonal_segment.h"
#include "../basic/match.h"

void init_cbs();

struct No_score_correction
{
	void operator()(int &score, int i, int query_anchor, int mult) const
	{}
};

struct Bias_correction : public std::vector<float>
{
	Bias_correction(const sequence &seq);
	void operator()(float &score, int i, int query_anchor, int mult) const
	{
		score += (*this)[query_anchor + i * mult];
	}
	int operator()(const Hsp &hsp) const;
	int operator()(const Diagonal_segment &d) const;
	std::vector<int8_t> int8;
};

#endif