#pragma once
#include "../../basic/value.h"

namespace Geo {

inline Loc i(Loc j, Loc d) {
	return d + j;
}

inline Loc j(Loc i, Loc d) {
	return i - d;
}

inline Loc diag_sub_matrix(Loc d, Loc i0, Loc j0) {
	return d + j0 - i0;
}

inline Loc rev_diag(Loc d, Loc qlen, Loc tlen) {
	return -d + qlen - tlen;
}

inline Loc min_diag(Loc qlen, Loc tlen) {
	return -(tlen - 1);
}

inline Loc max_diag(Loc qlen, Loc tlen) {
	return qlen - 1;
}

inline Loc clip_diag(Loc d, Loc qlen, Loc tlen) {
	return std::min(std::max(d, min_diag(qlen, tlen)), max_diag(qlen, tlen));
}

inline void assert_diag_bounds(Loc d, Loc qlen, Loc tlen) {
	assert(d >= min_diag(qlen, tlen) && d <= max_diag(qlen, tlen));
}

}