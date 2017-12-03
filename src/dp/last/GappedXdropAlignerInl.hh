// Copyright 2011, 2013, 2014 Martin C. Frith

// Modified by Benjamin Buchfink 2017/11/18

#ifndef GAPPED_XDROP_ALIGNER_INL_HH
#define GAPPED_XDROP_ALIGNER_INL_HH

#include <algorithm>
#include <cassert>
//#include <stdexcept>

namespace cbrc {

template<typename T, int N> T arrayMax(T (&array)[N]) {
  return *std::max_element(array, array + N);
}

template<typename T, int N> T arrayMin(T (&array)[N]) {
  return *std::min_element(array, array + N);
}

template<typename T> int maxIndex(T a, T b) {
  return b > a ? 1 : 0;
}

template<typename T> int maxIndex(T a, T b, T c) {
  return c > a ? maxIndex(b, c) + 1 : maxIndex(a, b);
}

template<typename T> int maxIndex(T a, T b, T c, T d) {
  return d > a ? maxIndex(b, c, d) + 1 : maxIndex(a, b, c);
}

template<typename T> int maxIndex(T a, T b, T c, T d, T e) {
  return e > a ? maxIndex(b, c, d, e) + 1 : maxIndex(a, b, c, d);
}

template<typename T> int maxIndex(T a, T b, T c, T d, T e, T f) {
  return f > a ? maxIndex(b, c, d, e, f) + 1 : maxIndex(a, b, c, d, e);
}

template<typename T> int maxIndex(T a, T b, T c, T d, T e, T f, T g) {
  return g > a ? maxIndex(b, c, d, e, f, g) + 1 : maxIndex(a, b, c, d, e, f);
}

template<typename T> T maxValue(T a, T b) {
  return std::max(a, b);
}

template<typename T> T maxValue(T a, T b, T c) {
  return maxValue(maxValue(a, b), c);
}

template<typename T>
T whichFrame(std::size_t antidiagonal, T frame0, T frame1, T frame2) {
  switch (antidiagonal % 3) {
    case 0: return frame1;  // the +1 frame
    case 1: return frame2;  // the -1 frame
    case 2: return frame0;
	default: return 0;
  }
}

inline bool isAffineGaps(int delExistenceCost, int delExtensionCost,
			 int insExistenceCost, int insExtensionCost,
			 int gapUnalignedCost) {
  return gapUnalignedCost >= delExtensionCost + insExtensionCost +
    std::max(delExistenceCost, insExistenceCost);
}

// The next two functions will stop the alignment at delimiters.  But
// this is not guaranteed if bestScore > INF / 2.  We could avoid this
// restriction by replacing -INF / 2 with bestScore - INF.

inline const int *finiteBeg(const int *beg, const int *end) {
  while (beg < end && *beg <= -INF / 2)
    ++beg;
  return beg;
}

inline const int *finiteEnd(const int *beg, const int *end) {
  while (end > beg && *(end-1) <= -INF / 2)
    --end;
  return end;
}

inline bool isDelimiter(uchar c, const int *scores) {
  return scores[c] <= -INF;
}

/*
inline void checkGappedXdropScore(int bestScore) {
  // If this happens, sentinels/delimiters might not work:
  if (bestScore > INF / 2)
    throw std::overflow_error("score got too high in gapped extension");
}
*/

inline void updateBest1(int &bestScore,
			std::size_t &bestAntidiagonal,
			std::size_t &bestSeq1position,
			int score,
			std::size_t antidiagonal,
			std::size_t seq1position) {
  if (score > bestScore) {
    bestScore = score;
    bestAntidiagonal = antidiagonal;
    bestSeq1position = seq1position;
  }
}

inline void GappedXdropAligner::updateBest(int &bestScore, int score,
                                           std::size_t antidiagonal,
                                           const int *x0, const int *x0base) {
  if (score > bestScore) {
    bestScore = score;
    bestAntidiagonal = antidiagonal;
    bestSeq1position = static_cast<std::size_t>(x0 - x0base);
  }
}

inline void updateMaxScoreDrop(int &maxScoreDrop,
                               std::size_t numCells, int maxMatchScore) {
  // If the current antidiagonal touches a sentinel/delimiter, then
  // maxMatches is the maximum possible number of matches starting
  // from the next antidiagonal.
  int maxMatches = static_cast<int>(numCells - 1);
  maxScoreDrop = std::min(maxScoreDrop, maxMatches * maxMatchScore - 1);
}

inline void updateFiniteEdges(std::size_t *maxSeq1begs,
                              std::size_t *minSeq1ends,
                              const int *x0base, const int *x0end,
                              std::size_t numCells) {
  const int *x0beg = x0end - numCells;

  maxSeq1begs[0] = maxSeq1begs[1] + 1;
  maxSeq1begs[1] = finiteBeg(x0beg, x0end) - x0base;

  minSeq1ends[0] = minSeq1ends[1];
  minSeq1ends[1] = finiteEnd(x0beg, x0end) - x0base + 1;
}

inline void updateFiniteEdges3(std::size_t *maxSeq1begs,
                               std::size_t *minSeq1ends,
                               const int *x0base, const int *x0end,
                               std::size_t numCells) {
  const int *x0beg = x0end - numCells;

  maxSeq1begs[0] = maxSeq1begs[1];
  maxSeq1begs[1] = maxSeq1begs[2];
  maxSeq1begs[2] = maxSeq1begs[3];
  maxSeq1begs[3] = maxSeq1begs[4] + 1;
  maxSeq1begs[4] = maxSeq1begs[5];
  maxSeq1begs[5] = maxSeq1begs[6];
  maxSeq1begs[6] = finiteBeg(x0beg, x0end) - x0base;

  minSeq1ends[0] = minSeq1ends[1];
  minSeq1ends[1] = minSeq1ends[2];
  minSeq1ends[2] = minSeq1ends[3];
  minSeq1ends[3] = minSeq1ends[4];
  minSeq1ends[4] = minSeq1ends[5] + 1;
  minSeq1ends[5] = minSeq1ends[6];
  minSeq1ends[6] = finiteEnd(x0beg, x0end) - x0base;
}

}

#endif
