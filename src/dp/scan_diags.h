#pragma once
#include "score_profile.h"

namespace DP {

void scan_diags128(const LongScoreProfile<int8_t>& qp, Sequence s, int d_begin, int j_begin, int j_end, int* out);
void scan_diags64(const LongScoreProfile<int8_t>& qp, Sequence s, int d_begin, int j_begin, int j_end, int* out);
void scan_diags(const LongScoreProfile<int8_t>& qp, Sequence s, int d_begin, int d_end, int j_begin, int j_end, int* out);
int diag_alignment(const int* s, int count);

}