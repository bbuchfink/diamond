// Copyright 2011, 2012 Martin C. Frith

// The algorithm is based on these recurrence formulas, for
// generalized affine gap costs.  For standard affine gap costs, set
// gup=infinity.
//
// gop = gapExistenceCost
// gep = gapExtensionCost
// gup = gapUnalignedCost
// F = frameshiftCost
//
// The 1st sequence: s(1), s(2), s(3), ...
// The  0 frame of the 2nd sequence: t(0, 1), t(0, 2), t(0, 3), ...
// The +1 frame of the 2nd sequence: t(1, 1), t(1, 2), t(1, 3), ...
// The -1 frame of the 2nd sequence: t(2, 1), t(2, 2), t(2, 3), ...
//
// frame(j)  =  (j+1) % 3
// index(j)  =  (j-1) / 3
// matchScore(i, j)  =  the score for aligning s(i) with t(frame(j), index(j)).
//
// Initialization:
// x(i, 0)  =  y(i, 0)  =  z(i, 0)  =  -INF  (for all i >= 0)
// x(i, 1)  =  y(i, 1)  =  z(i, 1)  =  -INF  (for all i >= 0)
// x(i, 2)  =  y(i, 2)  =  z(i, 2)  =  -INF  (for all i >= 0)
// x(i, 3)  =  y(i, 3)  =  z(i, 3)  =  -INF  (for all i >= 0)
// x(0, j)  =  y(0, j)  =  z(0, j)  =  -INF  (for all j >= 0)
// x(0, 2)  =  0
//
// Recurrence (i > 0 and j > 3):
// X(i, j)  =  max[ x(i-1, j-3), x(i-1, j-2) - F, x(i-1, j-4) - F ]
// Y(i, j)  =  max[ y(i-1, j) - gep, y(i-1, j-3) - gup ]
// Z(i, j)  =  max[ z(i, j-3) - gep, z(i-1, j-3) - gup ]
// b(i, j)  =  max[ X(i, j), Y(i, j), Z(i, j) ]
// x(i, j)  =  b(i, j) + matchScore(i, j)
// y(i, j)  =  max[ b(i, j) - gop, Y(i, j) ]
// z(i, j)  =  max[ b(i, j) - gop, Z(i, j) ]

// The recurrences are calculated antidiagonal-by-antidiagonal, where:
// antidiagonal  =  i*3 + j

// We store x(i, j), y(i, j), and z(i, j) in the following way.
// xScores: xxxxxoxxxxxxxxxx7xx8xx9xxAAxxBBxxCCxxDDDxxEEExxFFF...
// yScores: xxxxxxxxxxxxxxxx7xx8xx9xxAAxxBBxxCCxxDDDxxEEExxFFF...
// zScores: xxxxxxxxxxxxxxxx7xx8xx9xxAAxxBBxxCCxxDDDxxEEExxFFF...
// "o" indicates a cell with score = 0.
// "x" indicates a pad cell with score = -INF.
// "7", "8", "9", "A", etc. indicate cells in antidiagonal 7, 8, 9, 10, etc.
//
// We put 2 pad cells between antidiagonals.  This is sometimes
// necessary for forward frame-shifts, when we look-back by 7
// antidiagonals.

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"
#include <iostream>  // for debugging

namespace cbrc {

// Puts 7 "dummy" antidiagonals at the start, so that we can safely
// look-back from subsequent antidiagonals.
void GappedXdropAligner::init3() {
  scoreOrigins.resize(0);
  scoreEnds.resize(1);

  initAntidiagonal3(0, 0, 0);
  initAntidiagonal3(0, 2, 0);
  initAntidiagonal3(0, 4, 0);
  initAntidiagonal3(0, 6, 0);
  initAntidiagonal3(0, 8, 0);
  initAntidiagonal3(0, 10, 0);
  initAntidiagonal3(0, 12, 0);

  std::fill_n(xScores.begin(), 14, -INF);
  std::fill_n(yScores.begin(), 14, -INF);
  std::fill_n(zScores.begin(), 14, -INF);

  xScores[5] = 0;

  bestAntidiagonal = 8;
}

void GappedXdropAligner::initAntidiagonal3(std::size_t seq1beg,
                                           std::size_t scoreEnd,
                                           std::size_t numCells) {
  scoreOrigins.push_back(scoreEnd - seq1beg + 1);
  std::size_t newEnd = scoreEnd + numCells + 2;  // + 2 pad cells
  resizeScoresIfSmaller(newEnd);
  scoreEnds.push_back(newEnd);
}

// If seq2beg is the DNA coordinate relative to the start:
// seq1end = (antidiagonal - 8 - seq2beg) / 3 + 1
// seq2beg = antidiagonal - 8 - (seq1end - 1) * 3

// If the 0 frame is at the very end of the DNA sequence, then the +1
// frame will be just beyond a delimiter.  Which is OK.

// If the 0 frame is at the very start of the DNA sequence, then the
// -1 frame will be at an initial delimiter.  In that case, the code
// will miss alignments starting like this: reverse frameshift,
// deletion, insertion.  But it will find these equal-score
// alignments: reverse frameshift, insertion, deletion.

int GappedXdropAligner::align3(const uchar *seq1,
                               const uchar *seq2frame0,
                               const uchar *seq2frame1,  // the +1 frame
                               const uchar *seq2frame2,  // the -1 frame
                               bool isForward,
                               const ScoreMatrixRow *scorer,
                               int gapExistenceCost,
                               int gapExtensionCost,
                               int gapUnalignedCost,
                               int frameshiftCost,
                               int maxScoreDrop,
                               int maxMatchScore) {
  bool isAffine = gapUnalignedCost >= gapExistenceCost + 2 * gapExtensionCost;

  std::size_t maxSeq1begs[] = { 9, 9, 0, 9, 9, 9, 9 };
  std::size_t minSeq1ends[] = { 0, 0, 1, 0, 0, 0, 0 };

  int bestScore = 0;

  init3();

  for (std::size_t antidiagonal = 7; /* noop */; ++antidiagonal) {
    std::size_t seq1beg = arrayMin(maxSeq1begs);
    std::size_t seq1end = arrayMax(minSeq1ends);

    if (seq1beg >= seq1end) break;

    std::size_t scoreEnd = scoreEnds.back();
    std::size_t numCells = seq1end - seq1beg;

    initAntidiagonal3(seq1beg, scoreEnd, numCells);

    const uchar *seq2 =
        whichFrame(antidiagonal, seq2frame0, seq2frame1, seq2frame2);

    std::size_t seq2pos = (antidiagonal - 7) / 3 - seq1beg;

    const uchar *s1 = isForward ? seq1 + seq1beg : seq1 - seq1beg - 1;
    const uchar *s2 = isForward ? seq2 + seq2pos : seq2 - seq2pos - 1;

    if (isDelimiter(*s2, *scorer)) {
      // prevent forward frameshifts from jumping over delimiters:
      if (maxSeq1begs[1] == seq1beg) ++maxSeq1begs[1];
      // Update maxScoreDrop in some clever way?
      // But be careful if the -1 frame starts in an initial delimiter.
    }

    int minScore = bestScore - maxScoreDrop;

    int *x0 = &xScores[scoreEnd];
    int *y0 = &yScores[scoreEnd];
    int *z0 = &zScores[scoreEnd];
    const int *y3 = &yScores[hori3(antidiagonal, seq1beg)];
    const int *z3 = &zScores[vert3(antidiagonal, seq1beg)];
    const int *x6 = &xScores[diag3(antidiagonal, seq1beg)];
    const int *x5 = &xScores[diag3(antidiagonal + 1, seq1beg)];
    const int *x7 = &xScores[diag3(antidiagonal - 1, seq1beg)];

    *x0++ = *y0++ = *z0++ = -INF;  // add one pad cell

    const int *x0last = x0 + numCells;

    *x0++ = *y0++ = *z0++ = -INF;  // add one pad cell

    const int *x0base = x0 - seq1beg;

    if (isAffine) {
      if (isForward)
        while (1) {
          int s = maxValue(*x5, *x7);
          int x = maxValue(*x6, s - frameshiftCost);
          int y = *y3 - gapExtensionCost;
          int z = *z3 - gapExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
            updateBest(bestScore, b, antidiagonal, x0, x0base);
            *x0 = b + scorer[*s1][*s2];
            int g = b - gapExistenceCost;
            *y0 = maxValue(g, y);
            *z0 = maxValue(g, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          ++s1;  --s2;  ++x0;  ++y0;  ++z0;  ++y3;  ++z3;  ++x5;  ++x6;  ++x7;
        }
      else
        while (1) {
          int s = maxValue(*x5, *x7);
          int x = maxValue(*x6, s - frameshiftCost);
          int y = *y3 - gapExtensionCost;
          int z = *z3 - gapExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
			  
            updateBest(bestScore, b, antidiagonal, x0, x0base);
            *x0 = b + scorer[*s1][*s2];
            int g = b - gapExistenceCost;
            *y0 = maxValue(g, y);
            *z0 = maxValue(g, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          --s1;  ++s2;  ++x0;  ++y0;  ++z0;  ++y3;  ++z3;  ++x5;  ++x6;  ++x7;
        }
    } else {
      const int *y6 = &yScores[diag3(antidiagonal, seq1beg)];
      const int *z6 = &zScores[diag3(antidiagonal, seq1beg)];
      while (1) {
        int s = maxValue(*x5, *x7);
        int x = maxValue(*x6, s - frameshiftCost);
        int y = maxValue(*y3 - gapExtensionCost, *y6 - gapUnalignedCost);
        int z = maxValue(*z3 - gapExtensionCost, *z6 - gapUnalignedCost);
        int b = maxValue(x, y, z);
        if (b >= minScore) {
          updateBest(bestScore, b, antidiagonal, x0, x0base);
          *x0 = b + scorer[*s1][*s2];
          int g = b - gapExistenceCost;
          *y0 = maxValue(g, y);
          *z0 = maxValue(g, z);
        }
        else *x0 = *y0 = *z0 = -INF;
        if (x0 == x0last) break;
        ++x0;  ++y0;  ++z0;  ++y3;  ++z3;  ++x5;  ++x6;  ++x7;  ++y6;  ++z6;
        if (isForward) { ++s1;  --s2; }
        else           { --s1;  ++s2; }
      }
    }

    if (isDelimiter(*s1, *scorer))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    updateFiniteEdges3(maxSeq1begs, minSeq1ends, x0base, x0 + 1, numCells);
  }

  return bestScore;
}

bool GappedXdropAligner::getNextChunk3(std::size_t &end1,
                                       std::size_t &end2,
                                       std::size_t &length,
                                       int gapExistenceCost,
                                       int gapExtensionCost,
                                       int gapUnalignedCost,
                                       int frameshiftCost) {
  if (bestAntidiagonal == 8) return false;

  end1 = bestSeq1position;
  end2 = bestAntidiagonal - 8 - bestSeq1position * 3;
  length = 0;

  int state = 0;

  while (1) {
    if (state < 1 || state > 2) bestAntidiagonal -= 6;
    else                        bestAntidiagonal -= 3;

    if (state != 2) bestSeq1position -= 1;

    assert(bestAntidiagonal >= 7);
    assert(bestSeq1position * 3 <= bestAntidiagonal - 7);

    std::size_t h = hori3(bestAntidiagonal, bestSeq1position);
    std::size_t v = vert3(bestAntidiagonal, bestSeq1position);
    std::size_t d = diag3(bestAntidiagonal, bestSeq1position);
    std::size_t r = diag3(bestAntidiagonal + 1, bestSeq1position);
    std::size_t f = diag3(bestAntidiagonal - 1, bestSeq1position);

    int x = xScores[d];
    int y = yScores[h] - gapExtensionCost;
    int z = zScores[v] - gapExtensionCost;
    int a = yScores[d] - gapUnalignedCost;
    int b = zScores[d] - gapUnalignedCost;
    int i = xScores[r] - frameshiftCost;
    int j = xScores[f] - frameshiftCost;

    if (state == 1 || state == 5) {
      y += gapExistenceCost;
      a += gapExistenceCost;
    }

    if (state == 2 || state == 6) {
      z += gapExistenceCost;
      b += gapExistenceCost;
    }

    state = maxIndex(x, y, z, i, j, a, b);  // order?

    if (length == 0 && (state > 0 || bestAntidiagonal == 8))
      length = end1 - bestSeq1position;

    if (state == 3) {
      bestAntidiagonal += 1;
      state = 0;
    }

    if (state == 4) {
      bestAntidiagonal -= 1;
      state = 0;
    }

    if (length > 0 && state == 0) return true;
  }
}

}
