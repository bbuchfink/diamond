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

/** @file composition_adjustment.c
 * Highest level functions to solve the optimization problem for
 * compositional score matrix adjustment.
 *
 * @author Yi-Kuo Yu, Alejandro Schaffer, E. Michael Gertz
 */

#include <array>
#include <vector>
#include <algorithm>
#include <math.h>
#include "../basic/value.h"
#include "../basic/sequence.h"
#include "score_matrix.h"
#include "cbs.h"
#include "../basic/config.h"

namespace Stats {

using std::vector;

#define BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT    (1.e-5)
#define BLAST_KARLIN_LAMBDA_ITER_DEFAULT        17
#define COMPO_SCORE_MIN (-128)
#define LambdaRatioLowerBound 0.5
const int ALPH_TO_NCBI[] = { 1, 16, 13, 4, 3, 15, 5, 7, 8, 9, 11, 10, 12, 6, 14, 17, 18, 20, 22, 19 };

typedef struct Blast_ScoreFreq {
    int         score_min; /**< lowest allowed scores */
    int         score_max; /**< highest allowed scores */
    int         obs_min;   /**< lowest observed (actual) scores */
    int         obs_max;   /**< highest observed (actual) scores */
    double       score_avg; /**< average score, must be negative for local alignment. */
    double* sprob0;    /**< arrays for frequency of given score */
    double* sprob;     /**< arrays for frequency of given score, shifted down by score_min. */
} Blast_ScoreFreq;

inline double**
Nlm_DenseMatrixNew(int nrows,
    int ncols)
{
    int i;             /* iteration index */
    double** mat;     /* the new matrix */

    mat = (double**)calloc(nrows, sizeof(double*));
    if (mat != NULL) {
        mat[0] = (double*)malloc((size_t)nrows *
            (size_t)ncols * sizeof(double));
        if (mat[0] != NULL) {
            for (i = 1; i < nrows; i++) {
                mat[i] = &mat[0][i * ncols];
            }
        }
        else {
            free(mat);
            mat = NULL;
        }
    }
    return mat;
}

inline void
Nlm_DenseMatrixFree(double*** mat)
{
    if (*mat != NULL) {
        free((*mat)[0]);
        free(*mat);
    }
    *mat = NULL;
}

static int BLAST_Gcd(int a, int b)
{
    int   c;

    b = std::abs(b);
    if (b > a)
        c = a, a = b, b = c;

    while (b != 0) {
        c = a % b;
        a = b;
        b = c;
    }
    return a;
}

static double
NlmKarlinLambdaNR(double* probs, int d, int low, int high, double lambda0,
    double tolx, int itmax, int maxNewton, int* itn)
{
    int k;
    double x0, x, a = 0, b = 1;
    double f = 4;  /* Larger than any possible value of the poly in [0,1] */
    int isNewton = 0; /* we haven't yet taken a Newton step. */

    assert(d > 0);

    x0 = exp(-lambda0);
    x = (0 < x0 && x0 < 1) ? x0 : .5;

    for (k = 0; k < itmax; k++) { /* all iteration indices k */
        int i;
        double g, fold = f;
        int wasNewton = isNewton; /* If true, then the previous step was a */
                                  /* Newton step */
        isNewton = 0;            /* Assume that this step is not */

        /* Horner's rule for evaluating a polynomial and its derivative */
        g = 0;
        f = probs[low];
        for (i = low + d; i < 0; i += d) {
            g = x * g + f;
            f = f * x + probs[i];
        }
        g = x * g + f;
        f = f * x + probs[0] - 1;
        for (i = d; i <= high; i += d) {
            g = x * g + f;
            f = f * x + probs[i];
        }
        /* End Horner's rule */

        if (f > 0) {
            a = x; /* move the left endpoint */
        }
        else if (f < 0) {
            b = x; /* move the right endpoint */
        }
        else { /* f == 0 */
            break; /* x is an exact solution */
        }
        if (b - a < 2 * a * (1 - b) * tolx) {
            /* The midpoint of the interval converged */
            x = (a + b) / 2; break;
        }

        if (k >= maxNewton ||
            /* If convergence of Newton's method appears to be failing; or */
            (wasNewton && fabs(f) > .9 * fabs(fold)) ||
            /* if the previous iteration was a Newton step but didn't decrease
             * f sufficiently; or */
            g >= 0
            /* if a Newton step will move us away from the desired solution */
            ) { /* then */
          /* bisect */
            x = (a + b) / 2;
        }
        else {
            /* try a Newton step */
            double p = -f / g;
            double y = x + p;
            if (y <= a || y >= b) { /* The proposed iterate is not in (a,b) */
                x = (a + b) / 2;
            }
            else { /* The proposed iterate is in (a,b). Accept it. */
                isNewton = 1;
                x = y;
                if (fabs(p) < tolx * x * (1 - x)) break; /* Converged */
            } /* else the proposed iterate is in (a,b) */
        } /* else try a Newton step. */
    } /* end for all iteration indices k */
    *itn = k;
    return -log(x) / d;
}

double
Blast_KarlinLambdaNR(Blast_ScoreFreq* sfp, double initialLambdaGuess)
{
    int  low;        /* Lowest score (must be negative)  */
    int  high;       /* Highest score (must be positive) */
    int     itn;
    int i, d;
    double* sprob;
    double   returnValue;

    low = sfp->obs_min;
    high = sfp->obs_max;
    if (sfp->score_avg >= 0.) {   /* Expected score must be negative */
        return -1.0;
    }
    //if (BlastScoreChk(low, high) != 0) return -1.;

    sprob = sfp->sprob;
    /* Find greatest common divisor of all scores */
    for (i = 1, d = -low; i <= high - low && d > 1; ++i) {
        if (sprob[i + low] != 0.0) {
            d = BLAST_Gcd(d, i);
        }
    }
    returnValue =
        NlmKarlinLambdaNR(sprob, d, low, high,
            initialLambdaGuess,
            BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
            20, 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT, &itn);


    return returnValue;
}

double
s_CalcLambda(double probs[], int min_score, int max_score, double lambda0)
{

    int i;                 /* loop index */
    int score_range;       /* range of possible scores */
    double avg;            /* expected score of aligning two characters */
    Blast_ScoreFreq freq;  /* score frequency data */

    score_range = max_score - min_score + 1;
    avg = 0.0;
    for (i = 0; i < score_range; i++) {
        avg += (min_score + i) * probs[i];
    }
    freq.score_min = min_score;
    freq.score_max = max_score;
    freq.obs_min = min_score;
    freq.obs_max = max_score;
    freq.sprob0 = probs;
    freq.sprob = &probs[-min_score];
    freq.score_avg = avg;

    return Blast_KarlinLambdaNR(&freq, lambda0);
}

static void s_GetScoreRange(int* obs_min, int* obs_max,
    const int* const* matrix, int rows)
{
    int aa;                    /* index of an amino-acid in the 20
                                  letter alphabet */
    int irow, jcol;            /* matrix row and column indices */
    int minScore, maxScore;    /* largest and smallest observed scores */

    minScore = maxScore = 0;
    for (irow = 0; irow < rows; irow++) {
        for (aa = 0; aa < 20; aa++) {
            jcol = aa;
            if (matrix[irow][jcol] < minScore)
                minScore = matrix[irow][jcol];
            if (matrix[irow][jcol] > maxScore)
                maxScore = matrix[irow][jcol];
        }
    }
    *obs_min = minScore;
    *obs_max = maxScore;
}

int
s_GetMatrixScoreProbs(double** scoreProb, int* obs_min, int* obs_max,
    const int* const* matrix, int alphsize,
    const double* subjectProbArray,
    const double* queryProbArray)
{
    int aa;          /* index of an amino-acid in the 20 letter
                        alphabet */
    int irow, jcol;  /* matrix row and column indices */
    double* sprob;  /* a pointer to the element of the score
                        probabilities array that represents the
                        probability of the score 0*/
    int minScore;    /* smallest score in matrix; the same value as
                        (*obs_min). */
    int range;       /* the range of scores in the matrix */

    s_GetScoreRange(obs_min, obs_max, matrix, alphsize);
    minScore = *obs_min;
    range = *obs_max - *obs_min + 1;
    *scoreProb = (double*)calloc(range, sizeof(double));
    if (*scoreProb == NULL) {
        return -1;
    }
    sprob = &((*scoreProb)[-(*obs_min)]); /*center around 0*/
    for (irow = 0; irow < alphsize; irow++) {
        for (aa = 0; aa < 20; aa++) {
            jcol = aa;
            if (matrix[irow][jcol] >= minScore) {
                sprob[matrix[irow][jcol]] +=
                    (queryProbArray[irow] * subjectProbArray[jcol]);
            }
        }
    }
    return 0;
}

void
Blast_FreqRatioToScore(double** matrix, size_t rows, size_t cols, double Lambda)
{
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            if (0.0 == matrix[i][j]) {
                matrix[i][j] = COMPO_SCORE_MIN;
            }
            else {
            matrix[i][j] = log(matrix[i][j]) / Lambda;
            }
        }
    }
}

void
s_RoundScoreMatrix(int** matrix, size_t rows, size_t cols,
    double** floatScoreMatrix)
{
    for (size_t p = 0; p < rows; p++) {
        for (size_t c = 0; c < cols; c++) {
            if (floatScoreMatrix[p][c] < INT_MIN) {
                matrix[p][c] = INT_MIN;
            }
            else {
                matrix[p][c] = (int)std::round(floatScoreMatrix[p][c]);
            }
        }
    }
}


static double
s_CalcAvgScore(double* M, int alphsize, int incM, const double probs[])
{
    int j;                   /* iteration index */
    double score_iX = 0.0;   /* score of character i substituted by X */

    for (j = 0; j < alphsize; j++) {
        //if (alphaConvert[j] >= 0) {
            /* If the column corresponds to a true amino acid */
        score_iX += M[j * incM] * probs[j];
        //}
    }
    return score_iX;
}

static const double kMaximumXscore = -1.0;

static double
s_CalcXScore(double* M, int alphsize, int incM, const double probs[])
{
    return std::min(s_CalcAvgScore(M, alphsize, incM, probs), kMaximumXscore);
}

void
s_SetXUOScores(double** M, int alphsize,
    const double row_probs[], const double col_probs[])
{
    int i;                      /* iteration index */
    double score_XX = 0.0;      /* score of matching an X to an X */
    /* the matrix has alphsize colums (this variable exists just to
       make things easier to read) */
    const int cols = alphsize;

    for (i = 0; i < alphsize; i++) {
        //if (alphaConvert[i] >= 0) {
        double avg_iX = s_CalcAvgScore(M[i], alphsize, 1, col_probs);
        M[i][MASK_LETTER] = std::min(avg_iX, kMaximumXscore);
        score_XX += avg_iX * row_probs[i];

        M[MASK_LETTER][i] = s_CalcXScore(&M[0][i], alphsize, AMINO_ACID_COUNT, row_probs);
        //}
    }
    M[MASK_LETTER][MASK_LETTER] = std::min(score_XX, kMaximumXscore);

    /* Set X scores for pairwise ambiguity characters */
    //M[eBchar][eXchar] = s_CalcXScore(M[eBchar], alphsize, 1, col_probs);
    //M[eXchar][eBchar] = s_CalcXScore(&M[0][eBchar], alphsize, cols, row_probs);

    //M[eZchar][eXchar] = s_CalcXScore(M[eZchar], alphsize, 1, col_probs);
    //M[eXchar][eZchar] = s_CalcXScore(&M[0][eZchar], alphsize, cols, row_probs);
    //if (alphsize > eJchar) {
    //    M[eJchar][eXchar] = s_CalcXScore(M[eJchar], alphsize, 1, col_probs);
    //    M[eXchar][eJchar] =
    //        s_CalcXScore(&M[0][eJchar], alphsize, cols, row_probs);
    //}
    ///* Copy C scores to U */
    //memcpy(M[eSelenocysteine], M[eCchar], alphsize * sizeof(double));
    //for (i = 0; i < alphsize; i++) {
    //    M[i][eSelenocysteine] = M[i][eCchar];
    //}
    ///* Copy X scores to O */
    //if (alphsize > eOchar) {
    //    memcpy(M[eOchar], M[eXchar], alphsize * sizeof(double));
    //    for (i = 0; i < alphsize; i++) {
    //        M[i][eOchar] = M[i][eXchar];
    //    }
    //}
}

static int
s_ScaleSquareMatrix(int** matrix, int alphsize,
    const double row_prob[], const double col_prob[],
    double Lambda, const double(&freq_ratios)[NCBI_ALPH][NCBI_ALPH])
{
    double** scores;     /* a double precision matrix of scores */

    scores = Nlm_DenseMatrixNew(alphsize, alphsize);
    if (scores == 0) return -1;

    for (size_t i = 0; i < TRUE_AA; i++) {
        for (size_t j = 0; j < TRUE_AA; ++j)
            scores[i][j] = freq_ratios[ALPH_TO_NCBI[i]][ALPH_TO_NCBI[j]];
        //memcpy(scores[i], freq_ratios[i], alphsize * sizeof(double));
    }
    Blast_FreqRatioToScore(scores, TRUE_AA, TRUE_AA, Lambda);
    s_SetXUOScores(scores, TRUE_AA, row_prob, col_prob);
    s_RoundScoreMatrix(matrix, alphsize, alphsize, scores);
    /*for (i = 0; i < alphsize; i++) {
        matrix[i][(int)STOP_LETTER] = start_matrix[i][(int)STOP_LETTER];
        matrix[(int)STOP_LETTER][i] = start_matrix[(int)STOP_LETTER][i];
    }*/
    Nlm_DenseMatrixFree(&scores);

    return 0;
}

int
Blast_CompositionBasedStats(int** matrix, double* LambdaRatio,
    const int* const* matrix_in,
    const double queryProb[], const double resProb[], double lambda, const double(&freq_ratios)[NCBI_ALPH][NCBI_ALPH])
{
    double correctUngappedLambda; /* new value of ungapped lambda */
    int obs_min, obs_max;         /* smallest and largest score in the
                                     unscaled matrix */
    double* scoreArray;           /* an array of score probabilities */
    int out_of_memory;            /* status flag to indicate out of memory */

    out_of_memory = s_GetMatrixScoreProbs(&scoreArray, &obs_min, &obs_max, matrix_in, TRUE_AA, resProb, queryProb);
    const double ungappedLambda = lambda / config.cbs_matrix_scale;

    if (out_of_memory)
        return -1;
    correctUngappedLambda =
        s_CalcLambda(scoreArray, obs_min, obs_max, ungappedLambda);

    if (correctUngappedLambda < 0.0)
        return -1;

    /* calc_lambda will return -1 in the case where the
     * expected score is >=0; however, because of the MAX statement 3
     * lines below, LambdaRatio should always be > 0; the succeeding
     * test is retained as a vestige, in case one wishes to remove the
     * MAX statement and allow LambdaRatio to take on the error value
     * -1 */
    *LambdaRatio = correctUngappedLambda / ungappedLambda;
    //if (0 == pValueAdjustment)
    *LambdaRatio = std::min(1.0, *LambdaRatio);
    *LambdaRatio = std::max(*LambdaRatio, LambdaRatioLowerBound);

    if (*LambdaRatio > 0) {
        double scaledLambda = ungappedLambda / (*LambdaRatio);

        if (s_ScaleSquareMatrix(matrix, AMINO_ACID_COUNT,
            queryProb, resProb, scaledLambda, freq_ratios) < 0)
            return -1;

    }
    free(scoreArray);

    return 0;
}

vector<int> CompositionBasedStats(const int* const* matrix_in, const Composition& queryProb, const Composition& resProb, double lambda, const FreqRatios& freq_ratios) {
    vector<int> v(AMINO_ACID_COUNT*AMINO_ACID_COUNT);
    vector<int*> p;
    p.reserve(AMINO_ACID_COUNT);
    for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
        p.push_back(&v[i * AMINO_ACID_COUNT]);
    double LambdaRatio;
    if (Blast_CompositionBasedStats(p.data(), &LambdaRatio, matrix_in, queryProb.data(), resProb.data(), lambda, freq_ratios) != 0) {
        for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
            for (size_t j = 0; j < AMINO_ACID_COUNT; ++j)
                v[i * AMINO_ACID_COUNT + j] = score_matrix(i, j) * config.cbs_matrix_scale;
    }
    return v;
}

/* amino acid background frequencies from Robinson and Robinson */
static std::pair<char, double> Robinson_prob[] = {
      { 'A', 78.05 },
      { 'C', 19.25 },
      { 'D', 53.64 },
      { 'E', 62.95 },
      { 'F', 38.56 },
      { 'G', 73.77 },
      { 'H', 21.99 },
      { 'I', 51.42 },
      { 'K', 57.44 },
      { 'L', 90.19 },
      { 'M', 22.43 },
      { 'N', 44.87 },
      { 'P', 52.03 },
      { 'Q', 42.64 },
      { 'R', 51.29 },
      { 'S', 71.20 },
      { 'T', 58.41 },
      { 'V', 64.41 },
      { 'W', 13.30 },
      { 'Y', 32.16 }
}; /**< amino acid background frequencies from Robinson and Robinson */

double ideal_lambda(const int** matrix) {
    int obs_min, obs_max;
    double* scoreArray;
    double bg[TRUE_AA];
    double s = 0.0;
    for (size_t i = 0; i < TRUE_AA; ++i) {
        int j = value_traits.from_char(Robinson_prob[i].first);
        bg[j] = Robinson_prob[i].second;
        s += bg[j];
    }
    for (size_t i = 0; i < TRUE_AA; ++i)
        bg[i] /= s;

    int out_of_memory = s_GetMatrixScoreProbs(&scoreArray, &obs_min, &obs_max, matrix, TRUE_AA, bg, bg);

    if (out_of_memory)
        throw std::runtime_error("Failed lambda calculation.");
    double correctUngappedLambda = s_CalcLambda(scoreArray, obs_min, obs_max, 0.5);
    if (correctUngappedLambda < 0.0)
        throw std::runtime_error("Failed lambda calculation.");
    free(scoreArray);
    return correctUngappedLambda;
}

}