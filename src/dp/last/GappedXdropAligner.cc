// Copyright 2011, 2012, 2013 Martin C. Frith

// The algorithm is based on these recurrence formulas, for
// generalized affine gap costs.  For standard affine gap costs, set
// gup=infinity.
//
// gop = gapExistenceCost
// gep = gapExtensionCost
// gup = gapUnalignedCost
//
// The 1st sequence: s(1), s(2), s(3), ...
// The 2nd sequence: t(1), t(2), t(3), ...
//
// matchScore(i, j)  =  the score for aligning s(i) with t(j).
//
// Initialization:
// x(i, 0)  =  y(i, 0)  =  z(i, 0)  =  -INF  (for all i >= 0)
// x(0, j)  =  y(0, j)  =  z(0, j)  =  -INF  (for all j >= 0)
// x(0, 0)  =  0
//
// Recurrence (i > 0 and j > 0):
// X(i, j)  =  x(i-1, j-1)
// Y(i, j)  =  max[ y(i-1, j) - gep, y(i-1, j-1) - gup ]
// Z(i, j)  =  max[ z(i, j-1) - gep, z(i-1, j-1) - gup ]
// b(i, j)  =  max[ X(i, j), Y(i, j), Z(i, j) ]
// x(i, j)  =  b(i, j) + matchScore(i, j)
// y(i, j)  =  max[ b(i, j) - gop, Y(i, j) ]
// z(i, j)  =  max[ b(i, j) - gop, Z(i, j) ]
//
// Interpretations:
// X(i, j)  =  the best score for any alignment ending with s(i-1)
//             aligned to t(j-1).
// b(i, j)  =  the best score for any alignment ending at s(i-1) and
//             t(j-1).
// x(i, j)  =  the best score for any alignment ending with s(i)
//             aligned to t(j).

// The recurrences are calculated antidiagonal-by-antidiagonal, where:
// antidiagonal  =  i + j

// We store x(i, j), y(i, j), and z(i, j) in the following way.
// xScores: oxx2x33x444x5555x66666...
// yScores: xxx2x33x444x5555x66666...
// zScores: xxx2x33x444x5555x66666...
// "o" indicates a cell with score = 0.
// "x" indicates a pad cell with score = -INF.
// "2", "3", etc. indicate cells in antidiagonal 2, 3, etc.

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"
//#include <iostream>  // for debugging

namespace cbrc {

// Puts 2 "dummy" antidiagonals at the start, so that we can safely
// look-back from subsequent antidiagonals.
void GappedXdropAligner::init() {
  scoreOrigins.resize(0);
  scoreEnds.resize(1);

  initAntidiagonal(0, 0, 0);
  xScores[0] = 0;
  yScores[0] = -INF;
  zScores[0] = -INF;

  initAntidiagonal(0, 1, 0);
  xScores[1] = -INF;
  yScores[1] = -INF;
  zScores[1] = -INF;

  bestAntidiagonal = 0;
}

void GappedXdropAligner::initAntidiagonal(std::size_t seq1beg,
                                          std::size_t scoreEnd,
                                          std::size_t numCells) {
  scoreOrigins.push_back(scoreEnd - seq1beg);
  std::size_t newEnd = scoreEnd + numCells + 1;  // + 1 pad cell
  resizeScoresIfSmaller(newEnd);
  scoreEnds.push_back(newEnd);
}

int GappedXdropAligner::align(const uchar *seq1,
                              const uchar *seq2,
                              bool isForward,
			      int globality,
                              const ScoreMatrixRow *scorer,
                              int delExistenceCost,
                              int delExtensionCost,
                              int insExistenceCost,
                              int insExtensionCost,
                              int gapUnalignedCost,
                              int maxScoreDrop,
                              int maxMatchScore) {
  const bool isAffine = isAffineGaps(delExistenceCost, delExtensionCost,
				     insExistenceCost, insExtensionCost,
				     gapUnalignedCost);

  std::size_t maxSeq1begs[] = { 0, 9 };
  std::size_t minSeq1ends[] = { 1, 0 };

  int bestScore = 0;
  int bestEdgeScore = -INF;
  std::size_t bestEdgeAntidiagonal = 0;
  std::size_t bestEdgeSeq1position = 0;

  init();

  for (std::size_t antidiagonal = 0; /* noop */; ++antidiagonal) {
    std::size_t seq1beg = std::min(maxSeq1begs[0], maxSeq1begs[1]);
    std::size_t seq1end = std::max(minSeq1ends[0], minSeq1ends[1]);

    if (seq1beg >= seq1end) break;

    std::size_t scoreEnd = scoreEnds.back();
    std::size_t numCells = seq1end - seq1beg;

    initAntidiagonal(seq1beg, scoreEnd, numCells);

    std::size_t seq2pos = antidiagonal - seq1beg;

    const uchar *s1 = isForward ? seq1 + seq1beg : seq1 - seq1beg - 1;
    const uchar *s2 = isForward ? seq2 + seq2pos : seq2 - seq2pos - 1;

    if (!globality && isDelimiter(*s2, *scorer))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    int minScore = bestScore - maxScoreDrop;

    int *x0 = &xScores[scoreEnd];
    int *y0 = &yScores[scoreEnd];
    int *z0 = &zScores[scoreEnd];
    const int *y1 = &yScores[hori(antidiagonal, seq1beg)];
    const int *z1 = &zScores[vert(antidiagonal, seq1beg)];
    const int *x2 = &xScores[diag(antidiagonal, seq1beg)];

    const int *x0last = x0 + numCells;

    *x0++ = *y0++ = *z0++ = -INF;  // add one pad cell

    const int *x0base = x0 - seq1beg;

    if (globality && isDelimiter(*s2, *scorer)) {
      const int *z2 = &zScores[diag(antidiagonal, seq1beg)];
      int b = maxValue(*x2, *z1 - insExtensionCost, *z2 - gapUnalignedCost);
      if (b >= minScore)
	updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestEdgeSeq1position,
		    b, antidiagonal, seq1beg);
    }

    if (isAffine) {
      // We could avoid this code duplication, by using: transposed(scorer).
      if (isForward)
        while (1) {
          int x = *x2;
          int y = *y1 - delExtensionCost;
          int z = *z1 - insExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
            updateBest(bestScore, b, antidiagonal, x0, x0base);
            *x0 = b + scorer[*s1][*s2];
            *y0 = maxValue(b - delExistenceCost, y);
            *z0 = maxValue(b - insExistenceCost, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          ++s1;  --s2;  ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;
        }
      else
        while (1) {
          int x = *x2;
          int y = *y1 - delExtensionCost;
          int z = *z1 - insExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
            updateBest(bestScore, b, antidiagonal, x0, x0base);
            *x0 = b + scorer[*s1][*s2];
            *y0 = maxValue(b - delExistenceCost, y);
            *z0 = maxValue(b - insExistenceCost, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          --s1;  ++s2;  ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;
        }
    } else {
      const int *y2 = &yScores[diag(antidiagonal, seq1beg)];
      const int *z2 = &zScores[diag(antidiagonal, seq1beg)];
      while (1) {
        int x = *x2;
        int y = maxValue(*y1 - delExtensionCost, *y2 - gapUnalignedCost);
        int z = maxValue(*z1 - insExtensionCost, *z2 - gapUnalignedCost);
        int b = maxValue(x, y, z);
        if (b >= minScore) {
          updateBest(bestScore, b, antidiagonal, x0, x0base);
          *x0 = b + scorer[*s1][*s2];
          *y0 = maxValue(b - delExistenceCost, y);
          *z0 = maxValue(b - insExistenceCost, z);
        }
        else *x0 = *y0 = *z0 = -INF;
        if (x0 == x0last) break;
        ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;  ++y2;  ++z2;
        if (isForward) { ++s1;  --s2; }
        else           { --s1;  ++s2; }
      }
    }

    if (globality && isDelimiter(*s1, *scorer)) {
      const int *y2 = &yScores[diag(antidiagonal, seq1end-1)];
      int b = maxValue(*x2, *y1 - delExtensionCost, *y2 - gapUnalignedCost);
      if (b >= minScore)
	updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestEdgeSeq1position,
		    b, antidiagonal, seq1end-1);
    }

    if (!globality && isDelimiter(*s1, *scorer))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    updateFiniteEdges(maxSeq1begs, minSeq1ends, x0base, x0 + 1, numCells);
  }

  if (globality) {
    bestAntidiagonal = bestEdgeAntidiagonal;
    bestSeq1position = bestEdgeSeq1position;
    bestScore = bestEdgeScore;
  }
  return bestScore;
}

bool GappedXdropAligner::getNextChunk(std::size_t &end1,
                                      std::size_t &end2,
                                      std::size_t &length,
				      int delExistenceCost,
				      int delExtensionCost,
				      int insExistenceCost,
				      int insExtensionCost,
                                      int gapUnalignedCost) {
  if (bestAntidiagonal == 0) return false;

  end1 = bestSeq1position;
  end2 = bestAntidiagonal - bestSeq1position;
  const std::size_t undefined = -1;
  length = undefined;

  int state = 0;

  while (1) {
    assert(bestSeq1position <= bestAntidiagonal);

    std::size_t h = hori(bestAntidiagonal, bestSeq1position);
    std::size_t v = vert(bestAntidiagonal, bestSeq1position);
    std::size_t d = diag(bestAntidiagonal, bestSeq1position);

    int x = xScores[d];
    int y = yScores[h] - delExtensionCost;
    int z = zScores[v] - insExtensionCost;
    int a = yScores[d] - gapUnalignedCost;
    int b = zScores[d] - gapUnalignedCost;

    if (state == 1 || state == 3) {
      y += delExistenceCost;
      a += delExistenceCost;
    }

    if (state == 2 || state == 4) {
      z += insExistenceCost;
      b += insExistenceCost;
    }

    state = maxIndex(x, y, z, a, b);

    if (length == undefined && (state > 0 || bestAntidiagonal == 0)) {
      length = end1 - bestSeq1position;
      assert(length != undefined);
    }

    if (length != undefined && state == 0) return true;

    if (state < 1 || state > 2) bestAntidiagonal -= 2;
    else                        bestAntidiagonal -= 1;

    if (state != 2) bestSeq1position -= 1;
  }
}

}
