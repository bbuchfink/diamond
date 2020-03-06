#ifndef UNGAPPED_H_
#define UNGAPPED_H_

#include "../basic/value.h"
#include "../basic/diagonal_segment.h"
#include "comp_based_stats.h"
#include "../util/simd.h"

int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned seed_len, unsigned &delta, unsigned &len);
int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned &delta, unsigned &len);
int xdrop_ungapped_right(const Letter *query, const Letter *subject, int &len);
Diagonal_segment xdrop_ungapped(const sequence &query, const Bias_correction &query_bc, const sequence &subject, int qa, int sa);
Diagonal_segment xdrop_ungapped(const sequence &query, const sequence &subject, int qa, int sa);

namespace DP {

DECL_DISPATCH(void, window_ungapped, (const Letter *query, const Letter **subjects, size_t subject_count, int window, int *out))

}

#endif