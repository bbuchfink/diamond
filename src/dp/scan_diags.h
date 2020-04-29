#ifndef SCAN_DIAGS_H
#define SCAN_DIAGS_H

#include "score_profile.h"

namespace DP {

DECL_DISPATCH(void, scan_diags128, (const LongScoreProfile& qp, sequence s, int d_begin, int j_begin, int j_end, int* out))
DECL_DISPATCH(void, scan_diags64, (const LongScoreProfile& qp, sequence s, int d_begin, int j_begin, int j_end, int* out))
DECL_DISPATCH(int, diag_alignment, (const int* s, int count))

}

#endif