/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <array>
#include <vector>
#include "util/memory/memory_resource.h"
#include "basic/sequence.h"
#include "basic/statistics.h"
#include "standard_matrix.h"
#include "blast/matrix_adjust.h"

namespace Stats {

using Composition = std::array<MatrixFloat, TRUE_AA>;

Composition composition(const Sequence& s);

/** An collection of constants that specify all rules that may
 *  be used to generate a compositionally adjusted matrix.  */
typedef enum EMatrixAdjustRule {
    eDontAdjustMatrix = (-1),
    eCompoScaleOldMatrix = 0,
    eUnconstrainedRelEntropy = 1,
    eRelEntropyOldMatrixNewContext = 2,
    eRelEntropyOldMatrixOldContext = 3,
    eUserSpecifiedRelEntropy = 4
} EMatrixAdjustRule;

struct TargetMatrix {

    TargetMatrix(const Composition& query_comp, int query_len, unsigned cbs, const Sequence& target, Statistics& stats, std::pmr::monotonic_buffer_resource& pool, EMatrixAdjustRule rule);
    int score_width() const;

    std::pmr::vector<int8_t> scores;
    //std::vector<int32_t> scores32;
    //Scores<int8_t> scores_low, scores_high;

    int score_min, score_max;

};

/** Work arrays used to perform composition-based matrix adjustment */
typedef struct Blast_CompositionWorkspace {
    double** mat_b;       /**< joint probabilities for the matrix in
                                standard context */
    double** mat_final;   /**< optimized target frequencies */

    double* first_standard_freq;     /**< background frequency vector
                                           of the first sequence */
    double* second_standard_freq;    /**< background frequency vector of
                                           the second sequence */
} Blast_CompositionWorkspace;

/** Information about a amino-acid substitution matrix */
typedef struct Blast_MatrixInfo {
    char* matrixName;         /**< name of the matrix */
    int** startMatrix;     /**< Rescaled values of the original matrix */
    double** startFreqRatios;  /**< frequency ratios used to calculate matrix
                                    scores */
    int      rows;             /**< the number of rows in the scoring
                                    matrix. */
    int      cols;             /**< the number of columns in the scoring
                                    matrix, i.e. the alphabet size. */
    int      positionBased;    /**< is the matrix position-based */
    double   ungappedLambda;   /**< ungapped Lambda value for this matrix
                                    in standard context */
} Blast_MatrixInfo;

EMatrixAdjustRule adjust_matrix(const Composition& query_comp, int query_len, unsigned cbs, const Sequence& target);
void Blast_FreqRatioToScore(MatrixFloat** matrix, size_t rows, size_t cols, MatrixFloat Lambda);
void s_RoundScoreMatrix(int** matrix, size_t rows, size_t cols, MatrixFloat** floatScoreMatrix);
int s_GetMatrixScoreProbs(MatrixFloat** scoreProb, int* obs_min, int* obs_max,
    const int* const* matrix, int alphsize,
    const MatrixFloat* subjectProbArray,
    const MatrixFloat* queryProbArray);
MatrixFloat s_CalcLambda(MatrixFloat probs[], int min_score, int max_score, MatrixFloat lambda0);
MatrixFloat ideal_lambda(const int** matrix);
void s_SetXUOScores(MatrixFloat** M, int alphsize, const MatrixFloat row_probs[], const MatrixFloat col_probs[]);
int count_true_aa(const Sequence& s);
// bool use_seg_masking(const Sequence& a, const Sequence& b);

EMatrixAdjustRule
s_TestToApplyREAdjustmentConditional(int Len_query,
    int Len_match,
    const MatrixFloat* P_query,
    const MatrixFloat* P_match,
    const MatrixFloat* background_freqs);

struct CBS {
    static bool hauser(unsigned code) {
        switch (code) {
        case 0:
        case MATRIX_ADJUST:
        case CONDITIONAL_MATRIX_ADJUST:
        //case COMP_BASED_STATS:
        case COMP_BASED_STATS_AND_MATRIX_ADJUST:
        //case HAUSER_GLOBAL:
            return false;
        case 1:
        case 2:
        case HAUSER_AND_MATRIX_ADJUST:
            return true;
        default:
            throw std::runtime_error("Unknown CBS code.");
        }
    }
    static bool matrix_adjust(unsigned code) {
        switch (code) {
        case DISABLED:
        case HAUSER:
            return false;
        case DEPRECATED1:
        case HAUSER_AND_MATRIX_ADJUST:
        case MATRIX_ADJUST:
        case CONDITIONAL_MATRIX_ADJUST:
        //case COMP_BASED_STATS:
        case COMP_BASED_STATS_AND_MATRIX_ADJUST:
        //case HAUSER_GLOBAL:
            return true;
        default:
            throw std::runtime_error("Unknown CBS code.");
        }
    }
    static bool support_translated(unsigned code) {
        switch (code) {
        case DISABLED:
        case HAUSER:
            return true;
        default:
            return false;
        }
    }
    static bool conditioned(unsigned code) {
        switch (code) {
        case DEPRECATED1:
        case HAUSER_AND_MATRIX_ADJUST:
        case CONDITIONAL_MATRIX_ADJUST:
        case COMP_BASED_STATS_AND_MATRIX_ADJUST:
            return true;
        default:
            return false;
        }
    }
    static int tantan(unsigned code) {
        switch (code) {
        case DISABLED:
        case HAUSER:
            return 1;
        default:
            return 0;
        }
    }
    /*static int target_seg(unsigned code) {
        switch (code) {
        case DEPRECATED1:
        case HAUSER_AND_MATRIX_ADJUST:
        case MATRIX_ADJUST:
        case CONDITIONAL_MATRIX_ADJUST:
        //case COMP_BASED_STATS:
        //case COMP_BASED_STATS_AND_MATRIX_ADJUST:
        //case HAUSER_GLOBAL:
            return 1;
        default:
            return 0;
        }
    }*/
    enum {
        DISABLED = 0,
        HAUSER = 1,
        DEPRECATED1 = 2,
        HAUSER_AND_MATRIX_ADJUST = 3,
        MATRIX_ADJUST = 4,
		CONDITIONAL_MATRIX_ADJUST = 6,
        //COMP_BASED_STATS = 6,
        COMP_BASED_STATS_AND_MATRIX_ADJUST = 5,
        //HAUSER_GLOBAL = 8,
        COUNT
    };
    CBS(unsigned code, double query_match_distance_threshold, double length_ratio_threshold, double angle);
    double query_match_distance_threshold;
    double length_ratio_threshold;
    double angle;
};

void CompositionMatrixAdjust(int query_len, int target_len, const MatrixFloat* query_comp, const MatrixFloat* target_comp, int scale, MatrixFloat ungapped_lambda, const MatrixFloat* joint_probs, const MatrixFloat* background_freqs, std::array<int, AMINO_ACID_COUNT * AMINO_ACID_COUNT>& out, Statistics& stats);
bool CompositionBasedStats(const int* const* matrix_in, const Composition& queryProb, const Composition& resProb, double lambda, const FreqRatios& freq_ratios, std::array<int, AMINO_ACID_COUNT* AMINO_ACID_COUNT>& out);
std::vector<int> hauser_global(const Composition& query_comp, const Composition& target_comp);
int Blast_OptimizeTargetFrequencies(double x[],
    int alphsize,
    int* iterations,
    const double q[],
    const double row_sums[],
    const double col_sums[],
    int constrain_rel_entropy,
    double relative_entropy,
    double tol,
    int maxits);
bool OptimizeTargetFrequencies(double* out, const double* joints_prob, const double* row_probs, const double* col_probs, double relative_entropy, double tol, int maxits);

extern const int ALPH_TO_NCBI[];
extern CBS comp_based_stats;

}