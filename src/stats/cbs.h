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

#pragma once

#include <array>
#include <vector>
#include "../basic/sequence.h"
#include "standard_matrix.h"

namespace Stats {

using Composition = std::array<double, TRUE_AA>;

Composition composition(const Sequence& s);

struct TargetMatrix {

    TargetMatrix()
    {}

    TargetMatrix(const int16_t* query_matrix, const int16_t* target_matrix);

    TargetMatrix(const Composition& query_comp, int query_len, const Sequence& target);
    int score_width() const;
    bool blank() const {
        return scores.empty();
    }

    std::vector<int8_t> scores;
    std::vector<int32_t> scores32;

    int score_min, score_max;

};

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

void Blast_FreqRatioToScore(double** matrix, size_t rows, size_t cols, double Lambda);
void s_RoundScoreMatrix(int** matrix, size_t rows, size_t cols, double** floatScoreMatrix);
int s_GetMatrixScoreProbs(double** scoreProb, int* obs_min, int* obs_max,
    const int* const* matrix, int alphsize,
    const double* subjectProbArray,
    const double* queryProbArray);
double s_CalcLambda(double probs[], int min_score, int max_score, double lambda0);
double ideal_lambda(const int** matrix);
void s_SetXUOScores(double** M, int alphsize, const double row_probs[], const double col_probs[]);
int count_true_aa(const Sequence& s);
bool use_seg_masking(const Sequence& a, const Sequence& b);

EMatrixAdjustRule
s_TestToApplyREAdjustmentConditional(int Len_query,
    int Len_match,
    const double* P_query,
    const double* P_match,
    const double* background_freqs);

struct CBS {
    static bool hauser(unsigned code) {
        switch (code) {
        case 0:
        case 4:
        case COMP_BASED_STATS:
        case COMP_BASED_STATS_AND_MATRIX_ADJUST:
        case HAUSER_GLOBAL:
            return false;
        case 1:
        case 2:
        case 3:
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
        case COMP_BASED_STATS:
        case COMP_BASED_STATS_AND_MATRIX_ADJUST:
        case HAUSER_GLOBAL:
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
    static int target_seg(unsigned code) {
        switch (code) {
        case DEPRECATED1:
        case HAUSER_AND_MATRIX_ADJUST:
        case MATRIX_ADJUST:
        case COMP_BASED_STATS:
        case COMP_BASED_STATS_AND_MATRIX_ADJUST:
        case HAUSER_GLOBAL:
            return 1;
        default:
            return 0;
        }
    }
    enum {
        DISABLED = 0,
        HAUSER = 1,
        DEPRECATED1 = 2,
        HAUSER_AND_MATRIX_ADJUST = 3,
        MATRIX_ADJUST = 4,
        COMP_BASED_STATS = 5,
        COMP_BASED_STATS_AND_MATRIX_ADJUST = 6,
        HAUSER_GLOBAL = 7,
        COUNT
    };
    CBS(unsigned code, double query_match_distance_threshold, double length_ratio_threshold, double angle);
    double query_match_distance_threshold;
    double length_ratio_threshold;
    double angle;
    static constexpr int AVG_MATRIX_SCALE = 32;
};

std::vector<int> CompositionMatrixAdjust(int query_len, int target_len, const double* query_comp, const double* target_comp, int scale, double ungapped_lambda, const double* joint_probs, const double* background_freqs);
std::vector<int> CompositionBasedStats(const int* const* matrix_in, const Composition& queryProb, const Composition& resProb, double lambda, const FreqRatios& freq_ratios);
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

inline int16_t* make_16bit_matrix(const std::vector<int>& matrix) {
    int16_t* out = new int16_t[TRUE_AA * TRUE_AA];
    for (size_t i = 0; i < TRUE_AA; ++i)
        for (size_t j = 0; j < TRUE_AA; ++j)
            out[i * TRUE_AA + j] = matrix[i * AMINO_ACID_COUNT + j];
    return out;
}

extern const int ALPH_TO_NCBI[];
extern CBS comp_based_stats;

}