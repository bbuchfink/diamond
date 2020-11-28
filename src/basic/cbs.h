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
#include "sequence.h"

std::vector<double> composition(const sequence& s);

struct TargetMatrix {

    TargetMatrix()
    {}

    TargetMatrix(const int16_t* query_matrix, const int16_t* target_matrix);

    TargetMatrix(const double* query_comp, int query_len, const sequence& target);

    std::vector<int8_t> scores;
    std::vector<int32_t> scores32;

};

constexpr int TRUE_AA = 20;
extern const int ALPH_TO_NCBI[20];
extern const double BLOSUM62_FREQRATIOS[28][28];
constexpr double BLOSUM62_UNGAPPED_LAMBDA = 0.3176;
extern const double BLOSUM62_bg[20];

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

int
Blast_CompositionMatrixAdj(int** matrix,
    int alphsize,
    EMatrixAdjustRule matrix_adjust_rule,
    int length1,
    int length2,
    const double* stdaa_row_probs,
    const double* stdaa_col_probs,
    double lambda);

void Blast_FreqRatioToScore(double** matrix, int rows, int cols, double Lambda);
void s_RoundScoreMatrix(int** matrix, int rows, int cols, double** floatScoreMatrix);

EMatrixAdjustRule
s_TestToApplyREAdjustmentConditional(int Len_query,
    int Len_match,
    const double* P_query,
    const double* P_match);

struct CBS {
    static bool seg(unsigned code) {
        switch (code) {
        case 0:
        case 1:
            return false;
        case 2:
        case 3:
        case 4:
            return true;
        default:
            throw std::runtime_error("Unknown CBS code.");
        }
    }
    static bool hauser(unsigned code) {
        switch (code) {
        case 0:
        case 4:
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
        case HAUSER_AND_AVG_MATRIX_ADJUST:
        case HAUSER_AND_MATRIX_ADJUST:
        case MATRIX_ADJUST:
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
    static bool avg_matrix(unsigned code) {
        switch (code) {
        case HAUSER_AND_AVG_MATRIX_ADJUST:
            return true;
        default:
            return false;
        }
    }
    enum {
        DISABLED = 0,
        HAUSER = 1,
        HAUSER_AND_AVG_MATRIX_ADJUST = 2,
        HAUSER_AND_MATRIX_ADJUST = 3,
        MATRIX_ADJUST = 4,
        COUNT
    };
    static constexpr int AVG_MATRIX_SCALE = 32;
};

std::vector<int> CompositionMatrixAdjust(int query_len, int target_len, const double* query_comp, const double* target_comp, int scale);

inline int16_t* make_16bit_matrix(const std::vector<int>& matrix) {
    int16_t* out = new int16_t[TRUE_AA * TRUE_AA];
    for (int i = 0; i < TRUE_AA * TRUE_AA; ++i)
        out[i] = matrix[i];
    return out;
}