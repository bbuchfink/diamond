#pragma once

#include <array>
#include <vector>
#include "sequence.h"

std::array<double, 20> composition(const sequence& s);

struct TargetMatrix {

    TargetMatrix() :
        lambda_ratio(1.0)
    {}

    TargetMatrix(const double* query_comp, int query_len, const sequence& target);

    std::vector<int8_t> scores;
    std::vector<int32_t> scores32;
    double lambda_ratio;

};

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
    const double* stdaa_col_probs);

void Blast_FreqRatioToScore(double** matrix, int rows, int cols, double Lambda);
void s_RoundScoreMatrix(int** matrix, int rows, int cols, double** floatScoreMatrix);