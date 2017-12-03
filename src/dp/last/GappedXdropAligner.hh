// Copyright 2011, 2012, 2013 Martin C. Frith

// These routines extend an alignment in a given direction (forward or
// reverse) from given start points in two sequences.

// Forward alignments begin at the given start points, whereas reverse
// alignments begin one-before the given start points.

// To use: first call "align", which calculates the alignment but only
// returns its score.  To get the actual alignment, call
// "getNextChunk" to get each gapless chunk.

// The sequences had better end with sentinels: the score for matching
// a sentinel with anything should be -INF.

// The gap parameters correspond to "generalized affine gap costs"
// (Proteins 1998 32(1):88-96).
// gapExistenceCost = a
// gapExtensionCost = b
// gapUnalignedCost = c
// When c >= a + 2b, it reduces to standard affine gap costs:
// gap cost = gapExistenceCost + gapExtensionCost * (gap length).

// The insertion and deletion costs may differ.  Typically:
// delExistenceCost = insExistenceCost, and
// delExtensionCost = insExtensionCost.

// The algorithm proceeds antidiagonal-by-antidiagonal, similarly to
// Section 2 in J Comput Biol. 2000 7(1-2):203-14.  It does not allow
// the score to drop by more than maxScoreDrop below the highest score
// in any previous antidiagonal.

// If "globality" is 0, local alignment is performed: this finds the
// highest-scoring alignment ending anywhere.  Otherwise, overlap
// alignment is performed: this seeks the highest-scoring alignment
// ending at the end of either sequence.  If overlap alignment reaches
// the end of neither sequence (due to maxScoreDrop), -INF is
// returned.

// The parameter maxMatchScore should be the highest possible score
// for matching 2 letters.  This parameter is not actually necessary,
// but it provides some optimization opportunities.  If you give it a
// too-high value, the results will not change, but the run time may
// increase.

#ifndef GAPPED_XDROP_ALIGNER_HH
#define GAPPED_XDROP_ALIGNER_HH

#include "ScoreMatrixRow.hh"

#include <cstddef>  // size_t
#include <vector>

namespace cbrc {

typedef unsigned char uchar;

class TwoQualityScoreMatrix;

class GappedXdropAligner {
 public:
  int align(const uchar *seq1,  // start point in the 1st sequence
            const uchar *seq2,  // start point in the 2nd sequence
            bool isForward,  // forward or reverse extension?
	    int globality,
            const ScoreMatrixRow *scorer,  // the substitution score matrix
	    int delExistenceCost,
	    int delExtensionCost,
	    int insExistenceCost,
	    int insExtensionCost,
            int gapUnalignedCost,
            int maxScoreDrop,
            int maxMatchScore);

  // Like "align", but it aligns a sequence to a PSSM.
  int alignPssm(const uchar *seq,
                const ScoreMatrixRow *pssm,
                bool isForward,
		int globality,
		int delExistenceCost,
		int delExtensionCost,
		int insExistenceCost,
		int insExtensionCost,
                int gapUnalignedCost,
                int maxScoreDrop,
                int maxMatchScore);

  // Like "align", but both sequences have quality scores.
  int align2qual(const uchar *seq1,
                 const uchar *qual1,
                 const uchar *seq2,
                 const uchar *qual2,
                 bool isForward,
		 int globality,
                 const TwoQualityScoreMatrix &scorer,
		 int delExistenceCost,
		 int delExtensionCost,
		 int insExistenceCost,
		 int insExtensionCost,
                 int gapUnalignedCost,
                 int maxScoreDrop,
                 int maxMatchScore);

  // Call this repeatedly to get each gapless chunk of the alignment.
  // The chunks are returned in far-to-near order.  The chunk's end
  // coordinates in each sequence (relative to the start of extension)
  // and length are returned in the first 3 parameters.  If there are
  // no more chunks, the 3 parameters are unchanged and "false" is
  // returned.
  bool getNextChunk(std::size_t &end1,
                    std::size_t &end2,
                    std::size_t &length,
		    int delExistenceCost,
		    int delExtensionCost,
		    int insExistenceCost,
		    int insExtensionCost,
                    int gapUnalignedCost);

  // Like "align", but it aligns a protein sequence to a DNA sequence.
  // The DNA should be provided as 3 protein sequences, one for each
  // reading frame.  seq2frame0 is the in-frame start point.
  // seq2frame1 is shifted by 1 in the direction of alignment
  // extension.  seq2frame2 is shifted by 1 in the opposite direction.
  // The algorithm is "3-frame alignment", as described in J Comput
  // Biol. 1997 4(3):339-49.
  int align3(const uchar *seq1,
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
             int maxMatchScore);

  // Like "align3", but it aligns a protein sequence to a 3-frame PSSM.
  int align3pssm(const uchar *seq,
                 const ScoreMatrixRow *pssmFrame0,
                 const ScoreMatrixRow *pssmFrame1,
                 const ScoreMatrixRow *pssmFrame2,
                 bool isForward,
                 int gapExistenceCost,
                 int gapExtensionCost,
                 int gapUnalignedCost,
                 int frameshiftCost,
                 int maxScoreDrop,
                 int maxMatchScore);

  // Use this version of getNextChunk for protein-versus-DNA
  // alignment.  The end1 parameter receives a protein coordinate, and
  // end2 receives a DNA coordinate relative to the in-frame start.
  // The length parameter receives the number of residues/codons, not
  // bases.
  bool getNextChunk3(std::size_t &end1,
                     std::size_t &end2,
                     std::size_t &length,
                     int gapExistenceCost,
                     int gapExtensionCost,
                     int gapUnalignedCost,
                     int frameshiftCost);

  // The next 4 functions are for use by Centroid.  If the Centroid
  // code gets updated, it might make sense to change these functions too.

  // The number of antidiagonals, excluding dummy ones at the beginning.
  std::size_t numAntidiagonals() const
  { return scoreOrigins.size() - 2; }

  std::size_t scoreOrigin(std::size_t antidiagonal) const
  { return scoreOrigins[antidiagonal + 2]; }

  std::size_t numCellsAndPads(std::size_t antidiagonal) const
  { return scoreEnds[antidiagonal + 3] - scoreEnds[antidiagonal + 2]; }

  std::size_t scoreEndIndex(std::size_t antidiagonal) const
  { return scoreEnds[antidiagonal + 2]; }

  // The index in the score vectors, of the previous "horizontal" cell.
  std::size_t hori(std::size_t antidiagonal, std::size_t seq1coordinate) const
  { return scoreOrigins[antidiagonal + 1] + seq1coordinate; }

  // The index in the score vectors, of the previous "vertical" cell.
  std::size_t vert(std::size_t antidiagonal, std::size_t seq1coordinate) const
  { return scoreOrigins[antidiagonal + 1] + seq1coordinate + 1; }

  // The index in the score vectors, of the previous "diagonal" cell.
  std::size_t diag(std::size_t antidiagonal, std::size_t seq1coordinate) const
  { return scoreOrigins[antidiagonal] + seq1coordinate; }

  // The index in the score vectors, of the previous in-frame horizontal cell.
  std::size_t hori3(std::size_t antidiagonal, std::size_t seq1coordinate) const
  { return scoreOrigins[antidiagonal - 3] + seq1coordinate; }

  // The index in the score vectors, of the previous in-frame vertical cell.
  std::size_t vert3(std::size_t antidiagonal, std::size_t seq1coordinate) const
  { return scoreOrigins[antidiagonal - 3] + seq1coordinate + 1; }

  // The index in the score vectors, of the previous in-frame diagonal cell.
  std::size_t diag3(std::size_t antidiagonal, std::size_t seq1coordinate) const
  { return scoreOrigins[antidiagonal - 6] + seq1coordinate; }

  std::vector<int> xScores;  // best score ending with aligned letters
  std::vector<int> yScores;  // best score ending with insertion in seq1
  std::vector<int> zScores;  // best score ending with insertion in seq2

  std::vector<std::size_t> scoreOrigins;  // score origin for each antidiagonal
  std::vector<std::size_t> scoreEnds;  // score end pos for each antidiagonal

  // Our position during the trace-back:
  std::size_t bestAntidiagonal;
  std::size_t bestSeq1position;

  void resizeScoresIfSmaller(std::size_t size) {
    if (xScores.size() < size) {
      xScores.resize(size);
      yScores.resize(size);
      zScores.resize(size);
    }
  }

  void init();

  void initAntidiagonal(std::size_t seq1beg, std::size_t scoreEnd,
                        std::size_t numCells);

  void updateBest(int &bestScore, int score, std::size_t antidiagonal,
                  const int *x0, const int *x0ori);

  void init3();

  void initAntidiagonal3(std::size_t seq1beg, std::size_t scoreEnd,
                         std::size_t numCells);
};

}

#endif
