/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

/* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

#include "basic/config.h"
#include "score_matrix.h"
#include "blast/linear_algebra.h"
//#include "blast/nlm_linear_algebra.h"
//#include "linear_algebra.h"
#include "cbs.h"

using std::vector;
using std::array;

namespace Stats {

static constexpr int COMPO_NUM_TRUE_AA = 20;
static const int kReMatrixAdjustmentPseudocounts = 20;
/** relative entropy of BLOSUM62 */
static const double kFixedReBlosum62 = 0.44;

/* Documented in composition_adjustment.h. */
void
Blast_ApplyPseudocounts(MatrixFloat* probs20,
    int number_of_observations,
    const MatrixFloat* background_probs20)
{
    int i;                 /* loop index */
    MatrixFloat weight;         /* weight assigned to pseudocounts */
    MatrixFloat sum;            /* sum of the observed frequencies */
    /* pseudocounts as a double */
    MatrixFloat dpseudocounts = kReMatrixAdjustmentPseudocounts;
    /* Normalize probabilities */
    sum = 0.0;
    for (i = 0; i < COMPO_NUM_TRUE_AA; i++) {
        sum += probs20[i];
    }
    if (sum == 0.0) {  /* Can't normalize a zero vector */
        sum = 1.0;
    }
    weight = dpseudocounts / (number_of_observations + dpseudocounts);
    for (i = 0; i < COMPO_NUM_TRUE_AA; i++) {
        probs20[i] = (1.0 - weight) * probs20[i] / sum
            + weight * background_probs20[i];
    }
}

/* Documented in composition_adjustment.h. */
void
Blast_TrueAaToStdTargetFreqs(MatrixFloat** StdFreq, size_t StdAlphsize,
    const MatrixFloat* freq)
{
    /* Note I'm using a rough convention for this routine that uppercase
     * letters refer to quantities in the standard (larger) alphabet
     * and lowercase letters refer to the true amino acid (smaller)
     * alphabet.
     */
     /* Shorter names for the sizes of the two alphabets */
    const int small_alphsize = COMPO_NUM_TRUE_AA;
    size_t A, B;          /* characters in the std (big) alphabet */
    size_t a, b;          /* characters in the small alphabet */
    MatrixFloat sum;        /* sum of values in target_freq; used to normalize */
    sum = 0.0;
    for (a = 0; a < small_alphsize; a++) {
        for (b = 0; b < small_alphsize; b++) {
            sum += freq[a * TRUE_AA  + b];
        }
    }
    for (A = 0; A < StdAlphsize; A++) {
        /* for all rows */
        if (A >= TRUE_AA) {
            /* the row corresponds to a nonstandard reside */
            for (B = 0; B < StdAlphsize; B++) {
                StdFreq[A][B] = 0.0;
            }
        }
        else {
            /* the row corresponds to a standard reside */
            a = A;

            for (B = 0; B < StdAlphsize; B++) {
                /* for all columns */
                if (B >= TRUE_AA) {
                    /* the column corresponds to a nonstandard reside */
                    StdFreq[A][B] = 0.0;
                }
                else {
                    /* the column corresponds to a standard reside */
                    b = B;
                    StdFreq[A][B] = freq[a * TRUE_AA + b] / sum;
                }
            }
            /* Set values for two-character ambiguities */
            //StdFreq[A][eBchar] = StdFreq[A][eDchar] + StdFreq[A][eNchar];
            //StdFreq[A][eZchar] = StdFreq[A][eEchar] + StdFreq[A][eQchar];
            //if (StdAlphsize > eJchar) {
            //    StdFreq[A][eJchar] = StdFreq[A][eIchar] + StdFreq[A][eLchar];
            //}
        }
    }
    /* Add rows to set values for two-character ambiguities */
    //memcpy(StdFreq[eBchar], StdFreq[eDchar], StdAlphsize * sizeof(double));
    //Nlm_AddVectors(StdFreq[eBchar], StdAlphsize, 1.0, StdFreq[eNchar]);

    //memcpy(StdFreq[eZchar], StdFreq[eEchar], StdAlphsize * sizeof(double));
    //Nlm_AddVectors(StdFreq[eZchar], StdAlphsize, 1.0, StdFreq[eQchar]);

    //if (StdAlphsize > eJchar) {
    //    memcpy(StdFreq[eJchar], StdFreq[eIchar], StdAlphsize * sizeof(double));
    //    Nlm_AddVectors(StdFreq[eJchar], StdAlphsize, 1.0, StdFreq[eLchar]);
    //}
}

/* Documented in composition_adjustment.h. */
void
Blast_CalcFreqRatios(MatrixFloat** ratios, int alphsize,
    const MatrixFloat row_prob[], const MatrixFloat col_prob[])
{
    int i, j;
    for (i = 0; i < alphsize; i++) {
        if (row_prob[i] > 0) {
            for (j = 0; j < alphsize; j++) {
                if (col_prob[j] > 0) {
                    ratios[i][j] /= (row_prob[i] * col_prob[j]);
                }
            }
        }
    }
}

/**
 * Given a set of target frequencies and two sets of character
 * probabilities for the true amino acids in the ARND alphabet,
 * calculate a scoring matrix that has valid entries for all
 * characters in the NCBIstdaa amino acid alphabet.
 *
 * @param Matrix        the newly computed matrix
 * @param Alphsize      the size of the NCBIstdaa alphabet
 * @param target_freq   target frequencies for true amino acids (20x20)
 * @param StartMatrix   a matrix containing values for the stop character
 * @param row_prob      probabilities of true amino acids in the sequence
 *                      corresponding to the rows of matrix (length = 20)
 * @param col_prob      probabilities of true amino acids in the sequence
 *                      corresponding to the columns of matrix (length = 20)
 * @param Lambda        the desired scale of this matrix
 */
static int
s_ScoresStdAlphabet(int** Matrix, size_t Alphsize,
    const MatrixFloat* target_freq,
    const MatrixFloat row_prob[], const MatrixFloat col_prob[],
    MatrixFloat Lambda)
{
    /* Note: I'm using a rough convention for this routine that uppercase
     * letters refer to quantities in the standard (larger) alphabet
     * and lowercase letters refer to the true amino acid (smaller)
     * alphabet.
     */
    /* row and column probabilities in the NCBIstdaa alphabet */
    //double RowProb[COMPO_LARGEST_ALPHABET];
    //double ColProb[COMPO_LARGEST_ALPHABET];
    ///* A double precision score matrix */
    MatrixFloat** Scores = Nlm_DenseMatrixNew(Alphsize, Alphsize);
    if (Scores == NULL) {
        return -1;
    }
    //s_UnpackLetterProbs(RowProb, Alphsize, row_prob);
    //s_SetPairAmbigProbsToSum(RowProb, Alphsize);

    //s_UnpackLetterProbs(ColProb, Alphsize, col_prob);
    //s_SetPairAmbigProbsToSum(ColProb, Alphsize);

    Blast_TrueAaToStdTargetFreqs(Scores, Alphsize, target_freq);
    Blast_CalcFreqRatios(Scores, TRUE_AA, row_prob, col_prob);
    Blast_FreqRatioToScore(Scores, Alphsize, Alphsize, Lambda);
    s_SetXUOScores(Scores, TRUE_AA, row_prob, col_prob);

    s_RoundScoreMatrix(Matrix, Alphsize, Alphsize, Scores);
    Nlm_DenseMatrixFree(&Scores);

    //for (i = 0; i < Alphsize; i++) {
    //    Matrix[i][eStopChar] = StartMatrix[i][eStopChar];
    //    Matrix[eStopChar][i] = StartMatrix[eStopChar][i];
    //}
    return 0;
}

/* Documented in composition_adjustment.h. */
static int
Blast_CompositionMatrixAdj(int** matrix,
    EMatrixAdjustRule matrix_adjust_rule,
    int length1,
    int length2,
    const MatrixFloat* stdaa_row_probs,
    const MatrixFloat* stdaa_col_probs,
    double lambda,
    const MatrixFloat* joint_probs,
    const MatrixFloat* background_freqs,
    Statistics& stats)
{
    /*for (int i = 0; i < 20; ++i)
        printf("%.10f ", stdaa_row_probs[i]);
    printf("\n");
    for (int i = 0; i < 20; ++i)
        printf("%.10f ", stdaa_col_probs[i]);
    printf("\n");*/
    int iteration_count, status;
    MatrixFloat row_probs[COMPO_NUM_TRUE_AA], col_probs[COMPO_NUM_TRUE_AA];
    /* Target RE when optimizing the matrix; zero if the relative
       entropy should not be constrained. */
    MatrixFloat desired_re = 0.0;
    std::copy(stdaa_row_probs, stdaa_row_probs + 20, row_probs);
    std::copy(stdaa_col_probs, stdaa_col_probs + 20, col_probs);

    switch (matrix_adjust_rule) {
    //case eUnconstrainedRelEntropy:
    //    desired_re = 0.0;
    //    break;
    //case eRelEntropyOldMatrixNewContext:
    //    /* Calculate the desired re using the new marginal probs with
    //       the old matrix */
    //    status = Blast_EntropyOldFreqNewContext(&desired_re, &dummy,
    //        &iteration_count,
    //        NRrecord->mat_b,
    //        row_probs, col_probs);
    //    if (status < 0)     /* Error, e.g. memory */
    //        return status;
    //    else if (status > 0) /* we could not calculate the desired re */
    //        desired_re = 0.0; /* so, leave the re unconstrained */

    //    break;
    //case eRelEntropyOldMatrixOldContext:
    //    desired_re = Blast_TargetFreqEntropy(NRrecord->mat_b);
    //    break;
    case eUserSpecifiedRelEntropy:
        desired_re = kFixedReBlosum62;
        break;
    default:  /* I assert that we can't get here */
        fprintf(stderr, "Unknown flag for setting relative entropy"
            "in composition matrix adjustment");
        exit(1);
    }
    Blast_ApplyPseudocounts(row_probs, length1,
        background_freqs);
    Blast_ApplyPseudocounts(col_probs, length2,
        background_freqs);

    array<MatrixFloat, TRUE_AA * TRUE_AA> mat_final;
    status = Blast_OptimizeTargetFrequencies(mat_final.data(),
	COMPO_NUM_TRUE_AA,
        &iteration_count,
        joint_probs,
        row_probs, col_probs,
		(desired_re > 0.0),
        desired_re,
        config.cbs_err_tolerance,
        config.cbs_it_limit);

    if (status != 0)            /* Did not compute the target freqs */
        return status;

    return s_ScoresStdAlphabet(matrix, AMINO_ACID_COUNT, mat_final.data(), row_probs, col_probs, lambda);
}


/** 180 degrees in half a circle */
#define HALF_CIRCLE_DEGREES 180
/** some digits of PI */
#define PI 3.1415926543
/** @{ thresholds used to determine which composition mode to use */
#define HIGH_PAIR_THRESHOLD 0.4
#define LENGTH_LOWER_THRESHOLD 50
/** @} */


/** Return true if length > 50 and the two most frequent letters
 * occur a total of more that 40% of the time. */
static int
s_HighPairFrequencies(const MatrixFloat* letterProbs, int length)
{
    int i; /*index*/
    MatrixFloat max, second; /*two highest letter probabilities*/

    if (length <= LENGTH_LOWER_THRESHOLD) {
        return false;
    }
    max = 0;
    second = 0;
    for (i = 0; i < COMPO_NUM_TRUE_AA; i++) {
        if (letterProbs[i] > second) {
            second = letterProbs[i];
            if (letterProbs[i] > max) {
                second = max;
                max = letterProbs[i];
            }
        }
    }
    return (max + second) > HIGH_PAIR_THRESHOLD;
}

/**
 * Return true if either the query or the matching sequences
 * passes the test in s_HighPairFrequencies. */
static int
s_HighPairEitherSeq(const MatrixFloat* P_query, int length1,
    const MatrixFloat* P_match, int length2)
{
    int result1, result2;

    result1 = s_HighPairFrequencies(P_query, length1);
    result2 = s_HighPairFrequencies(P_match, length2);

    return result1 || result2;
}

/* Documented in composition_adjustment.h. */
MatrixFloat
Blast_GetRelativeEntropy(const MatrixFloat A[], const MatrixFloat B[])
{
    int i;                 /* loop index over letters */
    MatrixFloat temp;           /* intermediate term */
    MatrixFloat value = 0.0;    /* square of relative entropy */

    for (i = 0; i < COMPO_NUM_TRUE_AA; i++) {
        temp = (A[i] + B[i]) / 2;
        if (temp > 0) {
            if (A[i] > 0) {
                value += A[i] * log(A[i] / temp) / 2;
            }
            if (B[i] > 0) {
                value += B[i] * log(B[i] / temp) / 2;
            }
        }
    }
    if (value < 0) {             /* must be numerical rounding error */
        value = 0;
    }
    return sqrt(value);
}

/**
 * A function used to choose a mode for composition-based statistics.
 * Decide whether a relative-entropy score adjustment should be used
 * based on lengths and letter counts of the two matched sequences;
 * matrix_name is the underlying score matrix */
EMatrixAdjustRule
s_TestToApplyREAdjustmentConditional(int Len_query,
    int Len_match,
    const MatrixFloat* P_query,
    const MatrixFloat* P_match,
    const MatrixFloat* background_freqs)
{
    EMatrixAdjustRule which_rule; /* which relative entropy mode to
                                     return */
    int i;                       /* loop indices */
    MatrixFloat p_query[COMPO_NUM_TRUE_AA];
    MatrixFloat p_match[COMPO_NUM_TRUE_AA]; /*letter probabilities
                                                for query and match*/
    const MatrixFloat* p_matrix;       /* letter probabilities used in
                                     constructing matrix name*/
    double D_m_mat, D_q_mat, D_m_q;  /* distances between match and
                                        original between query and
                                        original between match and
                                        query*/
    double corr_factor = 0.0;     /* correlation between how p_query
                                     and p_match deviate from p_matrix
                                     */
    double len_q, len_m;          /* lengths of query and matching
                                     sequence in floating point */
    double len_large, len_small;  /* store the larger and smaller of
                                     len_q and len_m */
    double angle;                 /* angle between query and match
                                     probabilities */

    p_matrix = background_freqs; // Blast_GetMatrixBackgroundFreq(matrix_name);

    for (i = 0; i < COMPO_NUM_TRUE_AA; i++) {
        p_query[i] = P_query[i];
        p_match[i] = P_match[i];
        corr_factor +=
            (p_query[i] - p_matrix[i]) * (p_match[i] - p_matrix[i]);
    }
    D_m_mat = Blast_GetRelativeEntropy(p_match, p_matrix);
    D_q_mat = Blast_GetRelativeEntropy(p_query, p_matrix);
    D_m_q = Blast_GetRelativeEntropy(p_match, p_query);

    angle =
        acos((D_m_mat * D_m_mat + D_q_mat * D_q_mat -
            D_m_q * D_m_q) / 2.0 / D_m_mat / D_q_mat);
    /* convert from radians to degrees */
    angle = angle * HALF_CIRCLE_DEGREES / PI;

    len_q = 1.0 * Len_query;
    len_m = 1.0 * Len_match;
    if (len_q > len_m) {
        len_large = len_q;
        len_small = len_m;
    }
    else {
        len_large = len_m;
        len_small = len_q;
    }
    if (s_HighPairEitherSeq(P_query, Len_query, P_match, Len_match)) {
        which_rule = eUserSpecifiedRelEntropy;
    }
    else {
        if ((D_m_q > comp_based_stats.query_match_distance_threshold) &&
            (len_large / len_small > comp_based_stats.length_ratio_threshold) &&
            (angle > comp_based_stats.angle)) {
            which_rule = eCompoScaleOldMatrix;
        }
        else {
            which_rule = eUserSpecifiedRelEntropy;
        }
    }
    return which_rule;
}

void CompositionMatrixAdjust(int query_len, int target_len, const MatrixFloat* query_comp, const MatrixFloat* target_comp, int scale, MatrixFloat ungapped_lambda, const MatrixFloat* joint_probs, const MatrixFloat* background_freqs, array<int, AMINO_ACID_COUNT * AMINO_ACID_COUNT>& out, Statistics& stats) {
    array<int*, AMINO_ACID_COUNT> p;
    for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
        p[i] = &out[i * AMINO_ACID_COUNT];
    int r = Blast_CompositionMatrixAdj(p.data(),
        eUserSpecifiedRelEntropy,
        query_len,
        target_len,
        query_comp,
        target_comp,
        ungapped_lambda / scale,
        joint_probs,
        background_freqs, stats);
    if (r != 0) {
        for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
            for (size_t j = 0; j < AMINO_ACID_COUNT; ++j)
                out[i * AMINO_ACID_COUNT + j] = score_matrix(i, j) * scale;
        //throw std::runtime_error("Error computing composition matrix adjust.");
    }
}

}
