#include <iostream>
#include "blast_stat.h"
#include "blast_encoding.h"
#include "blast_def.h"
#include "blast_psi_priv.h"
#include "blast/ncbi_math.h"
#include "blast_setup.h"
#include "ncbi_std.h"
#include "blast_options.h"
#include "blast_encoding.h"



#define BLAST_SCORE_RANGE_MAX   (BLAST_SCORE_MAX - BLAST_SCORE_MIN) /**< maximum allowed range of BLAST scores. */
#define BLAST_KARLIN_LAMBDA0_DEFAULT    0.5 /**< Initial guess for the value of Lambda in BlastKarlinLambdaNR */
#define BLAST_KARLIN_K_SUMLIMIT_DEFAULT 0.0001 /**< K_SUMLIMIT_DEFAULT == sumlimit used in BlastKarlinLHtoK() */

#define BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT    (1.e-5) /**< LAMBDA_ACCURACY_DEFAULT == accuracy to which Lambda should be calc'd */

#define BLAST_KARLIN_LAMBDA_ITER_DEFAULT        17 /**< LAMBDA_ITER_DEFAULT == no. of iterations in LambdaBis = ln(accuracy)/ln(2)*/

#define BLAST_KARLIN_LAMBDA0_DEFAULT    0.5 /**< Initial guess for the value of Lambda in BlastKarlinLambdaNR */

#define BLAST_KARLIN_K_ITER_MAX 100 /**< upper limit on iterations for BlastKarlinLHtoK */

/** Number of statistical parameters in each row of the precomputed tables. */
#define BLAST_NUM_STAT_VALUES 11  /**< originally 8, now 11 to support Spouge's FSC. see notes below */



/** Holds values (gap-opening, extension, etc.) for a matrix. */
typedef double array_of_8[BLAST_NUM_STAT_VALUES];

/** Used to temporarily store matrix values for retrieval. */
typedef struct MatrixInfo {
    char*    name;       /**< name of matrix (e.g., BLOSUM90). */
    array_of_8  *values;    /**< The values (gap-opening, extension etc.). */
    Int4     *prefs;        /**< Preferences for display. */
    Int4     max_number_values;   /**< number of values (e.g., BLOSUM90_VALUES_MAX). */
} MatrixInfo;


/** Karlin-Altschul parameter values for substitution scores 1 and -5. */
static const array_of_8 blastn_values_1_5[] = {
        { 0, 0, 1.39, 0.747, 1.38, 1.00,  0, 100 },
        { 3, 3, 1.39, 0.747, 1.38, 1.00,  0, 100 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -4. */
static const array_of_8 blastn_values_1_4[] = {
        { 0, 0, 1.383, 0.738, 1.36, 1.02,  0, 100 },
        { 1, 2,  1.36,  0.67,  1.2,  1.1,  0,  98 },
        { 0, 2,  1.26,  0.43, 0.90,  1.4, -1,  91 },
        { 2, 1,  1.35,  0.61,  1.1,  1.2, -1,  98 },
        { 1, 1,  1.22,  0.35, 0.72,  1.7, -3,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -7.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
static const array_of_8 blastn_values_2_7[] = {
        { 0, 0,  0.69, 0.73, 1.34, 0.515,  0, 100 },
        { 2, 4,  0.68, 0.67,  1.2,  0.55,  0,  99 },
        { 0, 4,  0.63, 0.43, 0.90,   0.7, -1,  91 },
        { 4, 2, 0.675, 0.62,  1.1,   0.6, -1,  98 },
        { 2, 2,  0.61, 0.35, 0.72,   1.7, -3,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -3. */
static const array_of_8 blastn_values_1_3[] = {
        { 0, 0, 1.374, 0.711, 1.31, 1.05,  0, 100 },
        { 2, 2,  1.37,  0.70,  1.2,  1.1,  0,  99 },
        { 1, 2,  1.35,  0.64,  1.1,  1.2, -1,  98 },
        { 0, 2,  1.25,  0.42, 0.83,  1.5, -2,  91 },
        { 2, 1,  1.34,  0.60,  1.1,  1.2, -1,  97 },
        { 1, 1,  1.21,  0.34, 0.71,  1.7, -2,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -5.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
static const array_of_8 blastn_values_2_5[] = {
        { 0, 0, 0.675, 0.65,  1.1,  0.6, -1, 99 },
        { 2, 4,  0.67, 0.59,  1.1,  0.6, -1, 98 },
        { 0, 4,  0.62, 0.39, 0.78,  0.8, -2, 91 },
        { 4, 2,  0.67, 0.61,  1.0, 0.65, -2, 98 },
        { 2, 2,  0.56, 0.32, 0.59, 0.95, -4, 82 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -2. */
static const array_of_8 blastn_values_1_2[] = {
        { 0, 0, 1.28, 0.46, 0.85, 1.5, -2, 96 },
        { 2, 2, 1.33, 0.62,  1.1, 1.2,  0, 99 },
        { 1, 2, 1.30, 0.52, 0.93, 1.4, -2, 97 },
        { 0, 2, 1.19, 0.34, 0.66, 1.8, -3, 89 },
        { 3, 1, 1.32, 0.57,  1.0, 1.3, -1, 99 },
        { 2, 1, 1.29, 0.49, 0.92, 1.4, -1, 96 },
        { 1, 1, 1.14, 0.26, 0.52, 2.2, -5, 85 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -3.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
static const array_of_8 blastn_values_2_3[] = {
        { 0, 0,  0.55, 0.21, 0.46,  1.2, -5, 87 },
        { 4, 4,  0.63, 0.42, 0.84, 0.75, -2, 99 },
        { 2, 4, 0.615, 0.37, 0.72, 0.85, -3, 97 },
        { 0, 4,  0.55, 0.21, 0.46,  1.2, -5, 87 },
        { 3, 3, 0.615, 0.37, 0.68,  0.9, -3, 97 },
        { 6, 2,  0.63, 0.42, 0.84, 0.75, -2, 99 },
        { 5, 2, 0.625, 0.41, 0.78,  0.8, -2, 99 },
        { 4, 2,  0.61, 0.35, 0.68,  0.9, -3, 96 },
        { 2, 2, 0.515, 0.14, 0.33, 1.55, -9, 81 }
};

/** Karlin-Altschul parameter values for substitution scores 3 and -4. */
static const array_of_8 blastn_values_3_4[] = {
        { 6, 3, 0.389, 0.25, 0.56, 0.7, -5, 95},
        { 5, 3, 0.375, 0.21, 0.47, 0.8, -6, 92},
        { 4, 3, 0.351, 0.14, 0.35, 1.0, -9, 86},
        { 6, 2, 0.362, 0.16, 0.45, 0.8, -4, 88},
        { 5, 2, 0.330, 0.092, 0.28, 1.2, -13, 81},
        { 4, 2, 0.281, 0.046, 0.16, 1.8, -23, 69}
};

/** Karlin-Altschul parameter values for substitution scores 4 and -5. */
static const array_of_8 blastn_values_4_5[] = {
        { 0, 0, 0.22, 0.061, 0.22, 1.0, -15, 74 },
        { 6, 5, 0.28,  0.21, 0.47, 0.6 , -7, 93 },
        { 5, 5, 0.27,  0.17, 0.39, 0.7,  -9, 90 },
        { 4, 5, 0.25,  0.10, 0.31, 0.8, -10, 83 },
        { 3, 5, 0.23, 0.065, 0.25, 0.9, -11, 76 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -1. */
static const array_of_8 blastn_values_1_1[] = {
        { 3,  2, 1.09,  0.31, 0.55, 2.0,  -2, 99 },
        { 2,  2, 1.07,  0.27, 0.49, 2.2,  -3, 97 },
        { 1,  2, 1.02,  0.21, 0.36, 2.8,  -6, 92 },
        { 0,  2, 0.80, 0.064, 0.17, 4.8, -16, 72 },
        { 4,  1, 1.08,  0.28, 0.54, 2.0,  -2, 98 },
        { 3,  1, 1.06,  0.25, 0.46, 2.3,  -4, 96 },
        { 2,  1, 0.99,  0.17, 0.30, 3.3, -10, 90 }
};

/** Karlin-Altschul parameter values for substitution scores 3 and -2. */
static const array_of_8 blastn_values_3_2[] = {
        {  5,  5, 0.208, 0.030, 0.072, 2.9, -47, 77}
};

/** Karlin-Altschul parameter values for substitution scores 5 and -4. */
static const array_of_8 blastn_values_5_4[] = {
        { 10, 6, 0.163, 0.068, 0.16, 1.0, -19, 85 },
        {  8, 6, 0.146, 0.039, 0.11, 1.3, -29, 76 }
};



typedef struct BLAST_LetterProb {
    char  ch; /**< residue */
    double   p;  /**< probability of residue. */
} BLAST_LetterProb;

static BLAST_LetterProb nt_prob[] = {
        { 'A', 25.00 },
        { 'C', 25.00 },
        { 'G', 25.00 },
        { 'T', 25.00 }
}; /**< nucleotide probabilities (25% each letter) */


void**
_PSIDeallocateMatrix(void** matrix, unsigned int ncols)
{
    unsigned int i = 0;

    if (!matrix)
        return NULL;

    for (i = 0; i < ncols; i++) {
        sfree(matrix[i]);
    }
    sfree(matrix);
    return NULL;
}
static Int4
BlastKarlinEtoS_simple(double E, /* Expect value */
                       const Blast_KarlinBlk*  kbp,
                       Int8  searchsp)   /* size of search space */
{

    double   Lambda, K, H; /* parameters for Karlin statistics */
    Int4  S;
/* Smallest float that might not cause a floating point exception in
   S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda )); below.  */
    const double kSmallFloat = 1.0e-297;

    Lambda = kbp->Lambda;
    K = kbp->K;
    H = kbp->H;
    if (Lambda < 0. || K < 0. || H < 0.0)
    {
        return BLAST_SCORE_MIN;
    }

    E = MAX(E, kSmallFloat);

    S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda ));
    return S;
}

double
BLAST_KarlinStoE_simple(Int4 S,
                        Blast_KarlinBlk* kbp,
                        Int8  searchsp)   /* size of search space. */
{
    double   Lambda, K, H; /* parameters for Karlin statistics */

    Lambda = kbp->Lambda;
    K = kbp->K;
    H = kbp->H;
    if (Lambda < 0. || K < 0. || H < 0.) {
        return -1.;
    }

    return (double) searchsp * exp((double)(-Lambda * S) + kbp->logK);



}


static Blast_GumbelBlk*
s_BlastGumbelBlkNew() {
    return (Blast_GumbelBlk*) calloc(1, sizeof(Blast_GumbelBlk));
}
SBlastScoreMatrix*
SBlastScoreMatrixFree(SBlastScoreMatrix* matrix)
{
    if ( !matrix ) {
        return NULL;
    }

    if (matrix->data) {
        matrix->data = (int**) _PSIDeallocateMatrix((void**) matrix->data,
                                                    (unsigned)matrix->ncols);
    }

    /* Deallocate the matrix frequencies which is used by the
     * nucleotide custom matrix reader. -RMH-
     */
    if ( matrix->freqs )
        sfree(matrix->freqs);

    sfree(matrix);
    return NULL;
}
BlastScoringOptions*
BlastScoringOptionsFree(BlastScoringOptions* options)

{
    if (options == NULL)
        return NULL;

    sfree(options->matrix);
    sfree(options->matrix_path);
    sfree(options);

    return NULL;
}

Int4 BLAST_Gcd(Int4 a, Int4 b)
{
    Int4   c;

    b = ABS(b);
    if (b > a)
        c=a, a=b, b=c;

    while (b != 0) {
        c = a%b;
        a = b;
        b = c;
    }
    return a;
}

void**
_PSIAllocateMatrix(unsigned int ncols, unsigned int nrows,
                   unsigned int data_type_sz)
{
    void** retval = NULL;
    unsigned int i = 0;

    retval = (void**) malloc(sizeof(void*) * ncols);
    if ( !retval ) {
        return NULL;
    }

    for (i = 0; i < ncols; i++) {
        retval[i] = (void*) calloc(nrows, data_type_sz);
        if ( !retval[i] ) {
            retval = _PSIDeallocateMatrix(retval, i);
            break;
        }
    }
    return retval;
}


static Int2
BlastScoreBlkProteinMatrixRead(BlastScoreBlk* sbp, FILE *fp)
{
    char buf[512+3];
    char temp[512];
    char*   cp,* lp;
    char    ch;
    Int4 ** matrix;
    Int4 *  m;
    Int4 score;
    Uint4   a1cnt = 0, a2cnt = 0;
    char    a1chars[BLASTAA_SIZE], a2chars[BLASTAA_SIZE];
    long lineno = 0;
    double  xscore;
    register int  index1, index2;
    int x_index, u_index, o_index, c_index;
    const char kCommentChar = '#';
    const char* kTokenStr = " \t\n\r";

            ASSERT(sbp->alphabet_size == BLASTAA_SIZE);
            ASSERT(sbp->matrix);
            ASSERT(sbp->matrix->ncols == BLASTAA_SIZE);
            ASSERT(sbp->matrix->nrows == BLASTAA_SIZE);

    matrix = sbp->matrix->data;

    if (sbp->alphabet_code != BLASTNA_SEQ_CODE) {
        for (index1 = 0; index1 < sbp->alphabet_size; index1++)
            for (index2 = 0; index2 < sbp->alphabet_size; index2++)
                matrix[index1][index2] = BLAST_SCORE_MIN;
    }

    /* Read the residue names for the second alphabet */
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        ++lineno;
        if (strchr(buf, '\n') == NULL) {
            return 2;
        }

        if (buf[0] == kCommentChar) {
            /* save the comment line in a linked list */
            *strchr(buf, '\n') = NULLB;
            ListNodeCopyStr(&sbp->comments, 0, buf+1);
            continue;
        }
        if ((cp = strchr(buf, kCommentChar)) != NULL)
            *cp = NULLB;
        lp = (char*)strtok(buf, kTokenStr);
        if (lp == NULL) /* skip blank lines */
            continue;
        while (lp != NULL) {
            if (sbp->alphabet_code == BLASTAA_SEQ_CODE)
                ch = AMINOACID_TO_NCBISTDAA[toupper((unsigned char)(*lp))];
            else if (sbp->alphabet_code == BLASTNA_SEQ_CODE) {
                ch = IUPACNA_TO_BLASTNA[toupper((unsigned char)(*lp))];
            } else {
                ch = *lp;
            }
            a2chars[a2cnt++] = ch;
            lp = (char*)strtok(NULL, kTokenStr);
        }

        break; /* Exit loop after reading one line. */
    }

    if (a2cnt <= 1) {
        return 2;
    }

    while (fgets(buf, sizeof(buf), fp) != NULL)  {
        ++lineno;
        if ((cp = strchr(buf, '\n')) == NULL) {
            return 2;
        }
        if ((cp = strchr(buf, kCommentChar)) != NULL)
            *cp = NULLB;
        if ((lp = (char*)strtok(buf, kTokenStr)) == NULL)
            continue;
        ch = *lp;
        if ((cp = strtok(NULL, kTokenStr)) == NULL) {
            return 2;
        }
        if (a1cnt >= DIM(a1chars)) {
            return 2;
        }

        if (sbp->alphabet_code == BLASTAA_SEQ_CODE) {
            ch = AMINOACID_TO_NCBISTDAA[toupper((unsigned char) ch)];
        } else {
            if (sbp->alphabet_code == BLASTNA_SEQ_CODE) {
                ch = IUPACNA_TO_BLASTNA[toupper((unsigned char) ch)];
            }
        }
        a1chars[a1cnt++] = ch;
        m = &matrix[(int)ch][0];
        index2 = 0;
        while (cp != NULL) {
            if (index2 >= (int) a2cnt) {
                return 2;
            }
            strcpy(temp, cp);

            if (strcasecmp(temp, "na") == 0)  {
                score = BLAST_SCORE_MIN;
            } else  {
                if (sscanf(temp, "%lg", &xscore) != 1) {
                    return 2;
                }
                /*xscore = MAX(xscore, BLAST_SCORE_1MIN);*/
                if (xscore > BLAST_SCORE_MAX || xscore < BLAST_SCORE_MIN) {
                    return 2;
                }
                xscore += (xscore >= 0. ? 0.5 : -0.5);
                score = (Int4)xscore;
            }

            m[(int)a2chars[index2++]] = score;

            cp = strtok(NULL, kTokenStr);
        }
    }

    if (a1cnt <= 1) {
        return 2;
    }

    /* Use the C scores for U and X scores for O characters;
       if this is not done then they will never align to non-gap residues */
    x_index = AMINOACID_TO_NCBISTDAA['X'];
    u_index = AMINOACID_TO_NCBISTDAA['U'];
    o_index = AMINOACID_TO_NCBISTDAA['O'];
    c_index = AMINOACID_TO_NCBISTDAA['C'];
    for (index1 = 0; index1 < sbp->alphabet_size; index1++) {
        matrix[u_index][index1] = matrix[c_index][index1];
        matrix[index1][u_index] = matrix[index1][c_index];
        matrix[o_index][index1] = matrix[x_index][index1];
        matrix[index1][o_index] = matrix[index1][x_index];
    }

    return 0;
}


static Int2
BlastScoreBlkMaxScoreSet(BlastScoreBlk* sbp)
{
    Int4 score;
    Int4 ** matrix;
    Int2 index1, index2;

    sbp->loscore = BLAST_SCORE_MAX;
    sbp->hiscore = BLAST_SCORE_MIN;
    matrix = sbp->matrix->data;
    for (index1=0; index1<sbp->alphabet_size; index1++)
    {
        for (index2=0; index2<sbp->alphabet_size; index2++)
        {
            score = matrix[index1][index2];
            if (score <= BLAST_SCORE_MIN || score >= BLAST_SCORE_MAX)
                continue;
            if (sbp->loscore > score)
                sbp->loscore = score;
            if (sbp->hiscore < score)
                sbp->hiscore = score;
        }
    }
    /* If the lo/hi-scores are BLAST_SCORE_MIN/BLAST_SCORE_MAX, (i.e., for
    gaps), then use other scores. */

    if (sbp->loscore < BLAST_SCORE_MIN)
        sbp->loscore = BLAST_SCORE_MIN;
    if (sbp->hiscore > BLAST_SCORE_MAX)
        sbp->hiscore = BLAST_SCORE_MAX;

    return 0;
}


static Int2
BlastScoreBlkNucleotideMatrixRead(BlastScoreBlk* sbp, FILE *fp)
{
    Int4 ** matrix;
    double* freqs;
    Int4 i = 0;
    Int4 j = 0;
    Int4 rowIdx,colIdx,val,base;
    Int4 numFreqs = 0;
    Int4 alphaSize = 0;
    double fval;
    register int  index1, index2;
    char fbuf[512+3];
    char alphabet[24];
    char *cp,*ncp,*lp;
    double lambda_upper = 0;
    double lambda_lower = 0;
    double lambda = 0.5;
    double sum;
    double check;
    const char kCommentChar = '#';
    const char* kTokenStr = " \t\n\r";

    // Initialize matrix to default values
    matrix = sbp->matrix->data;
    for (index1 = 0; index1 < sbp->alphabet_size; index1++)
        for (index2 = 0; index2 < sbp->alphabet_size; index2++)
            matrix[index1][index2] = BLAST_SCORE_MIN;

    // Initialize matrix freqs to default values
    freqs = sbp->matrix->freqs;
    for (index1 = 0; index1 < sbp->alphabet_size; index1++)
        freqs[index1] = 0;

    alphabet[0] = 0;
    while ( fgets(fbuf,sizeof(fbuf),fp) ) {
        if (strchr(fbuf, '\n') == NULL) {
            return 2;
        }

        /* initialize column pointer */
        cp = (char *)fbuf;

        /* eat whitespace */
        while( (*cp) && isspace(*cp) ) cp++;

        if (*cp == kCommentChar) {
            /* special case FREQS line ( must exist ) */
            if ( (ncp = strstr( cp, (const char *)"FREQS" )) != NULL ) {
                cp = ncp + 5;
                /* eat whitespace */
                while( (*cp) && isspace(*cp) ) cp++;

                lp = (char*)strtok(cp, kTokenStr);
                /* Missing values */
                if (lp == NULL)
                    return 2;

                numFreqs = 0;
                while (lp != NULL) {
                    // Read Nucleotide
                    base = (int)IUPACNA_TO_BLASTNA[toupper((unsigned char)(*lp))];

                    lp = (char*)strtok(NULL, kTokenStr);
                    /* Expected a token pair */
                    if ( lp == NULL )
                        return 2;

                    // Read Frequency
                    if ( sscanf(lp, "%lf", &fval ) != 1 )
                        return( 2 );

                    // Store base/fval
                    freqs[base] = fval;
                    numFreqs++;
                    lp = (char*)strtok(NULL, kTokenStr);
                }
            }else {
                /* save the comment line in a linked list */
                *strchr(cp, '\n') = NULLB;
                ListNodeCopyStr(&sbp->comments, 0, cp);
            }
            continue;
        }

        /* alphabet line */
        if ( isalpha(*cp) && !alphabet[0] ) {
            j = 0;
            lp = (char*)strtok(cp, kTokenStr);
            while (lp != NULL) {
                alphabet[j++] = toupper((unsigned char)(*lp));
                lp = (char*)strtok(NULL, kTokenStr);
            }
            alphabet[j] = 0;
            alphaSize = j;
            continue;
        }else if ( isalpha(*cp) ) {
            /* Chew off first alphabet character */
            cp++;
            /* eat whitespace */
            while( (*cp) && isspace(*cp) ) cp++;
        }

        /* Matrix data */
        if ( isdigit(*cp) || *cp == '-' ) {
            j = 0;
            lp = (char*)strtok(cp, kTokenStr);
            rowIdx = (int)IUPACNA_TO_BLASTNA[toupper((unsigned char)alphabet[i])];
            while (lp != NULL) {
                if ( sscanf(lp, "%d", &val ) != 1 )
                    return( 2 );
                colIdx = (int)IUPACNA_TO_BLASTNA[toupper((unsigned char)alphabet[j++])];
                matrix[rowIdx][colIdx] = val;
                lp = (char*)strtok(NULL, kTokenStr);
            }
            /* We should have as many values as we do characters in the
               alphabet */
            if ( j != alphaSize )
                return( 2 );
            i++;
            continue;
        }
    }

    /* Expected 4 base frequencies, and a square matrix */
    if ( numFreqs != 4 || i != alphaSize )
        return( 2 );

    /* Calculate lambda for complexity adjusted scoring. This
       scoring system was designed by Phil Green and used in the
       cross_match package.  It was also used in MaskerAid to make
       wublast compatable with RepeatMasker. */
    do {
        sum = 0;
        check = 0;
        for ( i = 0 ; i < sbp->alphabet_size ; i++ ) {
            for ( j = 0 ; j < sbp->alphabet_size ; j++ ) {
                if ( freqs[i] && freqs[j] )
                {
                    sum += freqs[i] * freqs[j] *
                           exp( lambda * matrix[i][j] );
                    check += freqs[i] * freqs[j];
                }
            }
        }
                ASSERT( ( check < (double)1.001 ) && ( check > (double)0.999 ) );
        if ( sum < 1.0 ) {
            lambda_lower = lambda;
            lambda *= 2.0;
        }
    } while ( sum < 1.0 );

    lambda_upper = lambda;

    while ( lambda_upper - lambda_lower > (double).00001 ) {
        lambda = ( lambda_lower + lambda_upper ) / 2.0;
        sum          = 0;
        check        = 0;
        for ( i = 0 ; i < sbp->alphabet_size ; i++ )
        {
            for ( j = 0 ; j < sbp->alphabet_size ; j++ )
            {
                if ( freqs[i] && freqs[j] )
                {
                    sum += freqs[i] * freqs[j] *
                           exp( lambda * matrix[i][j] );
                    check += freqs[i] * freqs[j];
                }
            }
        }
                ASSERT( ( check < (double)1.001 ) && ( check > (double).999 ) );
        if ( sum >= 1.0 ) {
            lambda_upper = lambda;
        }
        else {
            lambda_lower = lambda;
        }
    }
    sbp->matrix->lambda = lambda;

    /* The value of 15 is a gap, which is a sentinel between strands in
       the ungapped extension algorithm. */
    for (index1=0; index1<BLASTNA_SIZE; index1++)
        matrix[BLASTNA_SIZE-1][index1] = INT4_MIN / 2;
    for (index1=0; index1<BLASTNA_SIZE; index1++)
        matrix[index1][BLASTNA_SIZE-1] = INT4_MIN / 2;

    return(0);

}

static const TNCBIScore s_IdentityPSM[25 * 25] = {
        /*       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,
                 F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *        */
        /*A*/    9, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*R*/   -5,  9, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*N*/   -5, -5,  9, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*D*/   -5, -5, -5,  9, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*C*/   -5, -5, -5, -5,  9, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*Q*/   -5, -5, -5, -5, -5,  9, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*E*/   -5, -5, -5, -5, -5, -5,  9, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*G*/   -5, -5, -5, -5, -5, -5, -5,  9, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*H*/   -5, -5, -5, -5, -5, -5, -5, -5,  9, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*I*/   -5, -5, -5, -5, -5, -5, -5, -5, -5,  9, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*L*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  9, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*K*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  9, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*M*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  9,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*F*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 9, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*P*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5,  9, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*S*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5,  9, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /*T*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5,  9, -5, -5, -5, -5, -5, -5, -5, -5,
        /*W*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5,  9, -5, -5, -5, -5, -5, -5, -5,
        /*Y*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5,  9, -5, -5, -5, -5, -5, -5,
        /*V*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5,  9, -5, -5, -5, -5, -5,
        /*B*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5,  9, -5, -5, -5, -5,
        /*J*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5,  9, -5, -5, -5,
        /*Z*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5,  9, -5, -5,
        /*X*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        /***/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
                 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  9
};
const SNCBIPackedScoreMatrix NCBISM_Identity = {
        "ARNDCQEGHILKMFPSTWYVBJZX*",
        s_IdentityPSM,
        -5
};

static const char kNCBIstdaa[] = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

int NCBISM_GetIndex(const SNCBIPackedScoreMatrix* sm, int aa)
{
    const char *p;

    /* Translate to NCBIeaa */
    if (aa >= 0  &&  aa < (int)sizeof(kNCBIstdaa)) {
        aa = kNCBIstdaa[aa];
    } else if (islower((unsigned char) aa)) {
        aa = toupper((unsigned char) aa);
    }

    p = strchr(sm->symbols, aa);
    return p ? (int)(p - sm->symbols) : -1;
}

TNCBIScore NCBISM_GetScore(const SNCBIPackedScoreMatrix* sm,
                           int aa1, int aa2)
{
    int i1, i2;
    i1 = NCBISM_GetIndex(sm, aa1);
    i2 = NCBISM_GetIndex(sm, aa2);
    if (i1 >=0  &&  i2 >= 0) {
        return sm->scores[(size_t)i1 * strlen(sm->symbols) + (size_t)i2];
    } else {
        return sm->defscore;
    }
}

long BLAST_Nint(double x)
{
    x += (x >= 0. ? 0.5 : -0.5);
    return (long)x;
}
Int2 BlastScoreBlkNuclMatrixCreate(BlastScoreBlk* sbp)
{
    Int2  index1, index2, degen;
    Int2 degeneracy[BLASTNA_SIZE+1];
    Int4 reward; /* reward for match of bases. */
    Int4 penalty; /* cost for mismatch of bases. */
    Int4** matrix; /* matrix to be populated. */
    /* How many of the first bases are ambiguous (four, of course). */
    const int k_number_non_ambig_bp = 4;

            ASSERT(sbp);
            ASSERT(sbp->alphabet_size == BLASTNA_SIZE);
            ASSERT(sbp->matrix);
            ASSERT(sbp->matrix->ncols == BLASTNA_SIZE);
            ASSERT(sbp->matrix->nrows == BLASTNA_SIZE);

    reward = sbp->reward;
    penalty = sbp->penalty;
    matrix = sbp->matrix->data;

    for (index1 = 0; index1<BLASTNA_SIZE; index1++)
        for (index2 = 0; index2<BLASTNA_SIZE; index2++)
            matrix[index1][index2] = 0;

    /* In blastna the 1st four bases are A, C, G, and T, exactly as it is
       ncbi2na. */
    /* ncbi4na gives them the value 1, 2, 4, and 8.  */
    /* Set the first four bases to degen. one */
    for (index1=0; index1<k_number_non_ambig_bp; index1++)
        degeneracy[index1] = 1;

    for (index1=k_number_non_ambig_bp; index1<BLASTNA_SIZE; index1++) {
        degen=0;
        for (index2=0; index2<k_number_non_ambig_bp; index2++) /* ncbi2na */
        {
            if (BLASTNA_TO_NCBI4NA[index1] & BLASTNA_TO_NCBI4NA[index2])
                degen++;
        }
        degeneracy[index1] = degen;
    }


    for (index1=0; index1<BLASTNA_SIZE; index1++) {
        for (index2=index1; index2<BLASTNA_SIZE; index2++) {
            if (BLASTNA_TO_NCBI4NA[index1] & BLASTNA_TO_NCBI4NA[index2]) {
                /* round up for positive scores, down for negatives. */
                matrix[index1][index2] =
                        BLAST_Nint( (double) ((degeneracy[index2]-1)*penalty +
                                              reward)/ (double) degeneracy[index2]);
                if (index1 != index2)
                {
                    matrix[index2][index1] = matrix[index1][index2];
                }
            }
            else
            {
                matrix[index1][index2] = penalty;
                matrix[index2][index1] = penalty;
            }
        }
    }

    /* The value of 15 is a gap, which is a sentinel between strands in
    the ungapped extension algorithm */
    for (index1=0; index1<BLASTNA_SIZE; index1++)
        matrix[BLASTNA_SIZE-1][index1] = INT4_MIN / 2;
    for (index1=0; index1<BLASTNA_SIZE; index1++)
        matrix[index1][BLASTNA_SIZE-1] = INT4_MIN / 2;

    return 0;
}

Int2 Blast_KarlinBlkCopy(Blast_KarlinBlk* kbp_to, Blast_KarlinBlk* kbp_from)
{
    if (!kbp_to || !kbp_from)
        return -1;

    kbp_to->Lambda = kbp_from->Lambda;
    kbp_to->K = kbp_from->K;
    kbp_to->logK = kbp_from->logK;
    kbp_to->H = kbp_from->H;
    kbp_to->paramC = kbp_from->paramC;
    return 0;
}
static
int /* bool */ s_NCBISM_StartsWith(const char* str, const char* pfx)
{
    for ( ;  *pfx;  ++str, ++pfx) {
        if (tolower((unsigned char)*str) != *pfx) {
            return 0;
        }
    }
    return 1;
}
const SNCBIPackedScoreMatrix* NCBISM_GetStandardMatrix(const char* name) {
    switch (name[0]) {
        case 'B':
        case 'b':
            if (!s_NCBISM_StartsWith(name, "blosum")) {
                return NULL;
            }


        case 'P':
        case 'p':
            if (!s_NCBISM_StartsWith(name, "pam")) {
                return NULL;
            }

        case 'I':
        case 'i':
            if (!s_NCBISM_StartsWith(name, "identity")) {
                return NULL;
            }
            return &NCBISM_Identity;

        default:
            return NULL;
    }
}

    static Int2
    BlastScoreBlkProteinMatrixLoad(BlastScoreBlk *sbp) {
        Int2 status = 0;
        Int4 **matrix = NULL;
        int i, j;   /* loop indices */
        int x_index, u_index, o_index, c_index;
        const SNCBIPackedScoreMatrix *psm;

                ASSERT(sbp);
        psm = NCBISM_GetStandardMatrix(sbp->name);
        if (psm == NULL)
            return 1;

                ASSERT(sbp->alphabet_size == BLASTAA_SIZE);
                ASSERT(sbp->matrix);
                ASSERT(sbp->matrix->ncols == BLASTAA_SIZE);
                ASSERT(sbp->matrix->nrows == BLASTAA_SIZE);

        matrix = sbp->matrix->data;

        /* Initialize with BLAST_SCORE_MIN */
        for (i = 0; i < sbp->alphabet_size; i++) {
            for (j = 0; j < sbp->alphabet_size; j++) {
                matrix[i][j] = BLAST_SCORE_MIN;
            }
        }

        for (i = 0; i < sbp->alphabet_size; i++) {
            for (j = 0; j < sbp->alphabet_size; j++) {
                /* skip special characters */
                if (i == AMINOACID_TO_NCBISTDAA['U'] ||
                    i == AMINOACID_TO_NCBISTDAA['O'] ||
                    i == AMINOACID_TO_NCBISTDAA['-'] ||
                    j == AMINOACID_TO_NCBISTDAA['U'] ||
                    j == AMINOACID_TO_NCBISTDAA['O'] ||
                    j == AMINOACID_TO_NCBISTDAA['-']) {
                    continue;
                }
                matrix[i][j] = NCBISM_GetScore((const SNCBIPackedScoreMatrix *) psm,
                                               i, j);
            }
        }

        /* Use the C scores for U and X scores for the O characters;
           if this is not done then they will never align to non-gap residues */
        x_index = AMINOACID_TO_NCBISTDAA['X'];
        u_index = AMINOACID_TO_NCBISTDAA['U'];
        o_index = AMINOACID_TO_NCBISTDAA['O'];
        c_index = AMINOACID_TO_NCBISTDAA['C'];
        for (i = 0; i < sbp->alphabet_size; i++) {
            matrix[u_index][i] = matrix[c_index][i];
            matrix[i][u_index] = matrix[i][c_index];
            matrix[o_index][i] = matrix[x_index][i];
            matrix[i][o_index] = matrix[i][x_index];
        }

        return status;
    }



Int2
Blast_ScoreBlkMatrixFill(BlastScoreBlk* sbp, GET_MATRIX_PATH get_path)
{
    bool matrix_found = FALSE;
    Int2 status = 0;

    /* For nucleotide case we first create a default matrix, based on the
       match and mismatch scores. */
    if (sbp->alphabet_code == BLASTNA_SEQ_CODE) {
        /* RMBLASTN supports reading a custom matrix.  Currently
         * I am using the sbp->read_in_matrix parameter to tell if
         * we should invoke the matrix reader.
         * -RMH-
         */
        if ( sbp->read_in_matrix && get_path )
        {
            // Read in custom rmblastn matrix
            matrix_found = FALSE;
        }
        else
        {
            if ( (status=BlastScoreBlkNuclMatrixCreate(sbp)) != 0)
                return status;
            matrix_found = TRUE;
        }
    }
    else
    {  /* Try to get compiled in matrix for proteins. */
        status = BlastScoreBlkProteinMatrixLoad(sbp);
        if (status == 0)
            matrix_found = TRUE;
    }

    /*if (matrix_found == FALSE && sbp->read_in_matrix && get_path) {
        // -RMH- Changed to FALSE
        char* matrix_path = get_path(sbp->name, FALSE);
        if (matrix_path) {

            FILE *fp = NULL;
            char* full_matrix_path = NULL;
            int path_len = strlen(matrix_path);
            int buflen = path_len + strlen(sbp->name);

            full_matrix_path = (char*) malloc((buflen + 1) * sizeof(char));
            if (!full_matrix_path) {
                return -1;
            }
            strncpy(full_matrix_path, matrix_path, buflen);
            strncat(full_matrix_path, sbp->name, buflen - path_len);

            sfree(matrix_path);

            if ( (fp=fopen(full_matrix_path, "r")) == NULL) {
                return -1;
            }
            sfree(full_matrix_path);

            if (sbp->alphabet_code == BLASTNA_SEQ_CODE) {
                // nucleotide blast ( rmblastn ) which can utilize a
                // a custom matrix -RMH-
                if ((status=BlastScoreBlkNucleotideMatrixRead(sbp, fp)) != 0)
                {
                    fclose(fp);
                    return status;
                }
            }else {
                // protein blast
                if ( (status=BlastScoreBlkProteinMatrixRead(sbp, fp)) != 0)
                {
                    fclose(fp);
                    return status;
                }
            }

            fclose(fp);
            matrix_found = TRUE;
        }
    }*/

    if (matrix_found == FALSE)
        return -1;

    if ( (status=BlastScoreBlkMaxScoreSet(sbp)) != 0)
        return status;

    return status;
}


SBlastScoreMatrix*
SBlastScoreMatrixNew(size_t ncols, size_t nrows)
{
    SBlastScoreMatrix* retval = NULL;

    retval = (SBlastScoreMatrix*) calloc(1, sizeof(SBlastScoreMatrix));
    if ( !retval ) {
        return SBlastScoreMatrixFree(retval);
    }

    retval->data = (int**) _PSIAllocateMatrix((unsigned)ncols, (unsigned)nrows, sizeof(int));
    if ( !retval->data ) {
        return SBlastScoreMatrixFree(retval);
    }

    /* Allocate additional attributes for use with custom
     * nucleotide matrices. -RMH-
     */
    retval->freqs = (double *) calloc(ncols, sizeof(double));
    retval->lambda = 0;

    retval->ncols = ncols;
    retval->nrows = nrows;
    return retval;
}

BlastScoreBlk*
BlastScoreBlkNew(Uint1 alphabet, Int4 number_of_contexts)

{
BlastScoreBlk* sbp;
char* use_old_fsc;

sbp = (BlastScoreBlk*) calloc(1, sizeof(BlastScoreBlk));

if ( !sbp ) {
return NULL;
}

sbp->alphabet_code = alphabet;
if (alphabet != BLASTNA_SEQ_CODE) {
sbp->alphabet_size = BLASTAA_SIZE;
} else {
sbp->alphabet_size = BLASTNA_SIZE;
}

/* set the alphabet type (protein or not). */
switch (alphabet) {
case BLASTAA_SEQ_CODE:
sbp->protein_alphabet = TRUE;
break;
case BLASTNA_SEQ_CODE:
sbp->protein_alphabet = FALSE;
break;
default:
break;
}

sbp->matrix = SBlastScoreMatrixNew(sbp->alphabet_size, sbp->alphabet_size);
if (sbp->matrix == NULL) {
return BlastScoreBlkFree(sbp);
}
sbp->scale_factor = 1.0;

/* FSCOLD: to switch back to original FSC, comment out the following line */
use_old_fsc = getenv("OLD_FSC");
if (!use_old_fsc) sbp->gbp = s_BlastGumbelBlkNew();

sbp->number_of_contexts = number_of_contexts;
sbp->sfp = (Blast_ScoreFreq**)
        calloc(sbp->number_of_contexts, sizeof(Blast_ScoreFreq*));
sbp->kbp_std = (Blast_KarlinBlk**)
        calloc(sbp->number_of_contexts, sizeof(Blast_KarlinBlk*));
sbp->kbp_gap_std = (Blast_KarlinBlk**)
        calloc(sbp->number_of_contexts, sizeof(Blast_KarlinBlk*));
sbp->kbp_psi = (Blast_KarlinBlk**)
        calloc(sbp->number_of_contexts, sizeof(Blast_KarlinBlk*));
sbp->kbp_gap_psi = (Blast_KarlinBlk**)
        calloc(sbp->number_of_contexts, sizeof(Blast_KarlinBlk*));

return sbp;
}
static Blast_GumbelBlk*
s_BlastGumbelBlkFree(Blast_GumbelBlk* gbp) {
    if ( !gbp) return NULL;
    sfree(gbp);
    return NULL;
}

Blast_KarlinBlk*
Blast_KarlinBlkFree(Blast_KarlinBlk* kbp)

{
    sfree(kbp);

    return kbp;
}

SPsiBlastScoreMatrix*
SPsiBlastScoreMatrixFree(SPsiBlastScoreMatrix* matrix)
{
    if ( !matrix ) {
        return NULL;
    }

    if (matrix->freq_ratios) {
        matrix->freq_ratios = (double**) _PSIDeallocateMatrix((void**)
                                                                      matrix->freq_ratios,
                                                              (unsigned)matrix->pssm->ncols);
    }

    matrix->pssm = SBlastScoreMatrixFree(matrix->pssm);
    matrix->kbp = Blast_KarlinBlkFree(matrix->kbp);
    sfree(matrix);
    return NULL;
}


Blast_ScoreFreq*
Blast_ScoreFreqFree(Blast_ScoreFreq* sfp)
{
    if (sfp == NULL)
        return NULL;

    if (sfp->sprob0 != NULL)
        sfree(sfp->sprob0);
    sfree(sfp);
    return sfp;
}
BlastScoreBlk*
BlastScoreBlkFree(BlastScoreBlk* sbp)

{
    Int4 index;
    if (sbp == NULL)
        return NULL;

    for (index=0; index<sbp->number_of_contexts; index++) {
        if (sbp->sfp)
            sbp->sfp[index] = Blast_ScoreFreqFree(sbp->sfp[index]);
        if (sbp->kbp_std)
            sbp->kbp_std[index] = Blast_KarlinBlkFree(sbp->kbp_std[index]);
        if (sbp->kbp_gap_std)
            sbp->kbp_gap_std[index] = Blast_KarlinBlkFree(sbp->kbp_gap_std[index]);
        if (sbp->kbp_psi)
            sbp->kbp_psi[index] = Blast_KarlinBlkFree(sbp->kbp_psi[index]);
        if (sbp->kbp_gap_psi)
            sbp->kbp_gap_psi[index] = Blast_KarlinBlkFree(sbp->kbp_gap_psi[index]);
    }
    if (sbp->kbp_ideal)
        sbp->kbp_ideal = Blast_KarlinBlkFree(sbp->kbp_ideal);
    if (sbp->gbp)
        sbp->gbp = s_BlastGumbelBlkFree(sbp->gbp);
    sfree(sbp->sfp);
    sbp->kbp = NULL;
    sbp->kbp_gap = NULL;
    sfree(sbp->kbp_std);
    sfree(sbp->kbp_psi);
    sfree(sbp->kbp_gap_std);
    sfree(sbp->kbp_gap_psi);
    sbp->matrix = SBlastScoreMatrixFree(sbp->matrix);
    sbp->comments = ListNodeFreeData(sbp->comments);
    if (sbp->name) {
        sfree(sbp->name);
    }
    sbp->psi_matrix = SPsiBlastScoreMatrixFree(sbp->psi_matrix);
    sfree(sbp->ambiguous_res);
    sfree(sbp);

    return NULL;
}
Blast_KarlinBlk*
Blast_KarlinBlkNew(void)

{
    Blast_KarlinBlk* kbp;

    kbp = (Blast_KarlinBlk*) calloc(1, sizeof(Blast_KarlinBlk));

    return kbp;
}
Blast_ResFreq*
Blast_ResFreqNew(const BlastScoreBlk* sbp)
{
    Blast_ResFreq* rfp;

    if (sbp == NULL)
    {
        return NULL;
    }

    rfp = (Blast_ResFreq*) calloc(1, sizeof(Blast_ResFreq));
    if (rfp == NULL)
        return NULL;

    rfp->alphabet_code = sbp->alphabet_code;

    rfp->prob0 = (double*) calloc(sbp->alphabet_size, sizeof(double));
    if (rfp->prob0 == NULL)
    {
        rfp = Blast_ResFreqFree(rfp);
        return rfp;
    }
    rfp->prob = rfp->prob0 - sbp->alphabet_start;

    return rfp;
}
static Int2
BlastScoreFreqCalc(const BlastScoreBlk* sbp, Blast_ScoreFreq* sfp, Blast_ResFreq* rfp1, Blast_ResFreq* rfp2)
{
    Int4 **  matrix;
    Int4  score, obs_min, obs_max;
    double      score_sum, score_avg;
    Int2     alphabet_start, alphabet_end, index1, index2;

    if (sbp == NULL || sfp == NULL)
        return 1;

    if (sbp->loscore < sfp->score_min || sbp->hiscore > sfp->score_max)
        return 1;

    for (score = sfp->score_min; score <= sfp->score_max; score++)
        sfp->sprob[score] = 0.0;

    matrix = sbp->matrix->data;

    alphabet_start = sbp->alphabet_start;
    alphabet_end = alphabet_start + sbp->alphabet_size;
    for (index1=alphabet_start; index1<alphabet_end; index1++)
    {
        for (index2=alphabet_start; index2<alphabet_end; index2++)
        {
            score = matrix[index1][index2];
            if (score >= sbp->loscore)
            {
                sfp->sprob[score] += rfp1->prob[index1] * rfp2->prob[index2];
            }
        }
    }

    score_sum = 0.;
    obs_min = obs_max = BLAST_SCORE_MIN;
    for (score = sfp->score_min; score <= sfp->score_max; score++)
    {
        if (sfp->sprob[score] > 0.)
        {
            score_sum += sfp->sprob[score];
            obs_max = score;
            if (obs_min == BLAST_SCORE_MIN)
                obs_min = score;
        }
    }
    sfp->obs_min = obs_min;
    sfp->obs_max = obs_max;

    score_avg = 0.0;
    if (score_sum > 0.0001 || score_sum < -0.0001)
    {
        for (score = obs_min; score <= obs_max; score++)
        {
            sfp->sprob[score] /= score_sum;
            score_avg += score * sfp->sprob[score];
        }
    }
    sfp->score_avg = score_avg;

    return 0;
}



Int2
Blast_ScoreBlkKbpIdealCalc(BlastScoreBlk* sbp)

{
    Blast_ResFreq* stdrfp = NULL;
    Blast_ScoreFreq* sfp = NULL;
    Int2 status = 0;

    if ( !sbp ) {
        return (status = 1);
    }

    stdrfp = Blast_ResFreqNew(sbp);
    Blast_ResFreqStdComp(sbp, stdrfp);
    sfp = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
    BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
    sbp->kbp_ideal = Blast_KarlinBlkNew();
    Blast_KarlinBlkUngappedCalc(sbp->kbp_ideal, sfp);

    stdrfp = Blast_ResFreqFree(stdrfp);
    sfp = Blast_ScoreFreqFree(sfp);

    return status;
}
Blast_ResFreq*
Blast_ResFreqFree(Blast_ResFreq* rfp)
{
    if (rfp == NULL)
        return NULL;

    if (rfp->prob0 != NULL)
        sfree(rfp->prob0);

    sfree(rfp);

    return rfp;
}

static double
BlastKarlinLHtoK(Blast_ScoreFreq* sfp, double lambda, double H)
{
    /*The next array stores the probabilities of getting each possible
      score in an alignment of fixed length; the array is shifted
      during part of the computation, so that
      entry 0 is for score 0.  */
    double         *alignmentScoreProbabilities = NULL;
    Int4            low;    /* Lowest score (must be negative) */
    Int4            high;   /* Highest score (must be positive) */
    Int4            range;  /* range of scores, computed as high - low*/
    double          K;      /* local copy of K  to return*/
    int             i;   /*loop index*/
    int             iterCounter; /*counter on iterations*/
    Int4            divisor; /*candidate divisor of all scores with
                               non-zero probabilities*/
    /*highest and lowest possible alignment scores for current length*/
    Int4            lowAlignmentScore, highAlignmentScore;
    Int4            first, last; /*loop indices for dynamic program*/
    register double innerSum;
    double          oldsum, oldsum2;  /* values of innerSum on previous
                                         iterations*/
    double          outerSum;        /* holds sum over j of (innerSum
                                        for iteration j/j)*/

    double          score_avg; /*average score*/
    /*first term to use in the closed form for the case where
      high == 1 or low == -1, but not both*/
    double          firstTermClosedForm;  /*usually store H/lambda*/
    int             iterlimit; /*upper limit on iterations*/
    double          sumlimit; /*lower limit on contributions
                                to sum over scores*/

    /*array of score probabilities reindexed so that low is at index 0*/
    double         *probArrayStartLow;

    /*pointers used in dynamic program*/
    double         *ptrP, *ptr1, *ptr2, *ptr1e;
    double          expMinusLambda; /*e^^(-Lambda) */

    if (lambda <= 0. || H <= 0.) {
        /* Theory dictates that H and lambda must be positive, so
         * return -1 to indicate an error */
        return -1.;
    }

    /*Karlin-Altschul theory works only if the expected score
      is negative*/
    if (sfp->score_avg >= 0.0) {
        return -1.;
    }

    low   = sfp->obs_min;
    high  = sfp->obs_max;
    range = high - low;

    probArrayStartLow = &sfp->sprob[low];
    /* Look for the greatest common divisor ("delta" in Appendix of PNAS 87 of
       Karlin&Altschul (1990) */
    for (i = 1, divisor = -low; i <= range && divisor > 1; ++i) {
        if (probArrayStartLow[i] != 0.0)
            divisor = BLAST_Gcd(divisor, i);
    }

    high   /= divisor;
    low    /= divisor;
    lambda *= divisor;

    range = high - low;

    firstTermClosedForm = H/lambda;
    expMinusLambda      = exp((double) -lambda);

    if (low == -1 && high == 1) {
        K = (sfp->sprob[low*divisor] - sfp->sprob[high*divisor]) *
            (sfp->sprob[low*divisor] - sfp->sprob[high*divisor]) / sfp->sprob[low*divisor];
        return(K);
    }

    if (low == -1 || high == 1) {
        if (high != 1) {
            score_avg = sfp->score_avg / divisor;
            firstTermClosedForm
                    = (score_avg * score_avg) / firstTermClosedForm;
        }
        return firstTermClosedForm * (1.0 - expMinusLambda);
    }

    sumlimit  = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;
    iterlimit = BLAST_KARLIN_K_ITER_MAX;

    alignmentScoreProbabilities =
            (double *)calloc((iterlimit*range + 1), sizeof(*alignmentScoreProbabilities));
    if (alignmentScoreProbabilities == NULL)
        return -1.;

    outerSum = 0.;
    lowAlignmentScore = highAlignmentScore = 0;
    alignmentScoreProbabilities[0] = innerSum = oldsum = oldsum2 = 1.;

    for (iterCounter = 0;
         ((iterCounter < iterlimit) && (innerSum > sumlimit));
         outerSum += innerSum /= ++iterCounter) {
        first = last = range;
        lowAlignmentScore  += low;
        highAlignmentScore += high;
        /*dynamic program to compute P(i,j)*/
        for (ptrP = alignmentScoreProbabilities +
                    (highAlignmentScore-lowAlignmentScore);
             ptrP >= alignmentScoreProbabilities;
             *ptrP-- =innerSum) {
            ptr1  = ptrP - first;
            ptr1e = ptrP - last;
            ptr2  = probArrayStartLow + first;
            for (innerSum = 0.; ptr1 >= ptr1e; ) {
                innerSum += *ptr1  *  *ptr2;
                ptr1--;
                ptr2++;
            }
            if (first)
                --first;
            if (ptrP - alignmentScoreProbabilities <= range)
                --last;
        }
        /* Horner's rule */
        innerSum = *++ptrP;
        for( i = lowAlignmentScore + 1; i < 0; i++ ) {
            innerSum = *++ptrP + innerSum * expMinusLambda;
        }
        innerSum *= expMinusLambda;

        for (; i <= highAlignmentScore; ++i)
            innerSum += *++ptrP;
        oldsum2 = oldsum;
        oldsum  = innerSum;
    }




#ifdef ADD_GEOMETRIC_TERMS_TO_K
    /*old code assumed that the later terms in sum were
      asymptotically comparable to those of a geometric
      progression, and tried to speed up convergence by
      guessing the estimated ratio between sucessive terms
      and using the explicit terms of a geometric progression
      to speed up convergence. However, the assumption does not
      always hold, and convergenece of the above code is fast
      enough in practice*/
    /* Terms of geometric progression added for correction */
    {
        double     ratio;  /* fraction used to generate the
                                   geometric progression */

        ratio = oldsum / oldsum2;
        if (ratio >= (1.0 - sumlimit*0.001)) {
            K = -1.;
            if (alignmentScoreProbabilities != NULL)
                sfree(alignmentScoreProbabilities);
            return K;
        }
        sumlimit *= 0.01;
        while (innerSum > sumlimit) {
            oldsum   *= ratio;
            outerSum += innerSum = oldsum / ++iterCounter;
        }
    }
#endif

    K = -exp((double)-2.0*outerSum) /
        (firstTermClosedForm*BLAST_Expm1(-(double)lambda));

    if (alignmentScoreProbabilities != NULL)
        sfree(alignmentScoreProbabilities);

    return K;
}


static Int2
BlastScoreChk(Int4 lo, Int4 hi)
{
    if (lo >= 0 || hi <= 0 ||
        lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX)
        return 1;

    if (hi - lo > BLAST_SCORE_RANGE_MAX)
        return 1;

    return 0;
}

double BLAST_Expm1(double	x)
{
    double	absx = ABS(x);

    if (absx > .33)
        return exp(x) - 1.;

    if (absx < 1.e-16)
        return x;

    return x * (1. + x *
                     (1./2. + x *
                              (1./6. + x *
                                       (1./24. + x *
                                                 (1./120. + x *
                                                            (1./720. + x *
                                                                       (1./5040. + x *
                                                                                   (1./40320. + x *
                                                                                                (1./362880. + x *
                                                                                                              (1./3628800. + x *
                                                                                                                             (1./39916800. + x *
                                                                                                                                             (1./479001600. +
                                                                                                                                              x/6227020800.))))))))))));
}

double BLAST_Powi(double x, Int4 n)
{
    double   y;

    if (n == 0)
        return 1.;

    if (x == 0.) {
        if (n < 0) {
            return HUGE_VAL;
        }
        return 0.;
    }

    if (n < 0) {
        x = 1./x;
        n = -n;
    }

    y = 1.;
    while (n > 0) {
        if (n & 1)
            y *= x;
        n /= 2;
        x *= x;
    }
    return y;
}

static Int2
s_AdjustGapParametersByGcd(array_of_8* normal, array_of_8* linear, int size, Int4* gap_existence_max, Int4* gap_extend_max, int divisor)
{
    if (divisor == 1)
        return 0;

    if (size <=0)
        return 1;

    (*gap_existence_max) *= divisor;
    (*gap_extend_max) *= divisor;

    if (normal)
    {
        int i;

        for (i=0; i<size; i++)
        {  /* divide lambda and alpha by divisor. */
            /* multiply gap existence and extension by divisor. */
            normal[i][0] *= divisor;
            normal[i][1] *= divisor;
            normal[i][2] /= divisor;
            normal[i][5] /= divisor;
        }
    }
    if (linear)
    {  /* divide lambda and alpha by divisor. */
        linear[0][0] *= divisor;
        linear[0][1] *= divisor;
        linear[0][2] /= divisor;
        linear[0][5] /= divisor;
    }

    return 0;
}


static Int2
s_SplitArrayOf8(const array_of_8* input, const array_of_8** normal, const array_of_8** non_affine, bool *split)
{

    if (input == NULL || normal == NULL || non_affine == NULL)
        return -1;

    *normal = NULL;
    *non_affine = NULL;

    if (input[0][0] == 0 && input[0][1] == 0)
    {
        *normal = input+1;
        *non_affine = input;
        *split = TRUE;
    }
    else
    {
        *normal = input;
        *split = FALSE;
    }
    return 0;

}


static Int2
s_GetNuclValuesArray(Int4 reward, Int4 penalty, Int4* array_size,
                     array_of_8** normal, array_of_8** non_affine,
                     Int4* gap_open_max, Int4* gap_extend_max, bool* round_down,
                     Blast_Message** error_return)
{
    Int2 status = 0;
    const array_of_8 * kValues = NULL;
    const array_of_8 * kValues_non_affine = NULL;
    bool split = FALSE;
    int divisor = BLAST_Gcd(reward, penalty);

    *round_down = FALSE;

    *array_size = 0;
    *normal = NULL;
    *non_affine = NULL;

    if (divisor != 1)
    {
        reward /= divisor;
        penalty /= divisor;
    }

    if (reward == 1 && penalty == -5) {
        if ((status=s_SplitArrayOf8(blastn_values_1_5, &kValues, &kValues_non_affine, &split)))
            return status;

        *array_size = sizeof(blastn_values_1_5)/sizeof(array_of_8);
        *gap_open_max = 3;
        *gap_extend_max = 3;
    } else if (reward == 1 && penalty == -4) {
        if ((status=s_SplitArrayOf8(blastn_values_1_4, &kValues, &kValues_non_affine, &split)))
            return status;

        *array_size = sizeof(blastn_values_1_4)/sizeof(array_of_8);
        *gap_open_max = 2;
        *gap_extend_max = 2;
    } else if (reward == 2 && penalty == -7) {
        if ((status=s_SplitArrayOf8(blastn_values_2_7, &kValues, &kValues_non_affine, &split)))
            return status;

        *round_down = TRUE;
        *array_size = sizeof(blastn_values_2_7)/sizeof(array_of_8);
        *gap_open_max = 4;
        *gap_extend_max = 4;
    } else if (reward == 1 && penalty == -3) {
        if ((status=s_SplitArrayOf8(blastn_values_1_3, &kValues, &kValues_non_affine, &split)))
            return status;

        *array_size = sizeof(blastn_values_1_3)/sizeof(array_of_8);
        *gap_open_max = 2;
        *gap_extend_max = 2;
    } else if (reward == 2 && penalty == -5) {
        if ((status=s_SplitArrayOf8(blastn_values_2_5, &kValues, &kValues_non_affine, &split)))
            return status;

        *round_down = TRUE;
        *array_size = sizeof(blastn_values_2_5)/sizeof(array_of_8);
        *gap_open_max = 4;
        *gap_extend_max = 4;
    } else if (reward == 1 && penalty == -2) {
        if ((status=s_SplitArrayOf8(blastn_values_1_2, &kValues, &kValues_non_affine, &split)))
            return status;

        *array_size = sizeof(blastn_values_1_2)/sizeof(array_of_8);
        *gap_open_max = 2;
        *gap_extend_max = 2;
    } else if (reward == 2 && penalty == -3) {
        if ((status=s_SplitArrayOf8(blastn_values_2_3, &kValues, &kValues_non_affine, &split)))
            return status;

        *round_down = TRUE;
        *array_size = sizeof(blastn_values_2_3)/sizeof(array_of_8);
        *gap_open_max = 6;
        *gap_extend_max = 4;
    } else if (reward == 3 && penalty == -4) {
        if ((status=s_SplitArrayOf8(blastn_values_3_4, &kValues, &kValues_non_affine, &split)))
            return status;

        *round_down = TRUE;
        *array_size = sizeof(blastn_values_3_4)/sizeof(array_of_8);
        *gap_open_max = 6;
        *gap_extend_max = 3;
    } else if (reward == 1 && penalty == -1) {
        if ((status=s_SplitArrayOf8(blastn_values_1_1, &kValues, &kValues_non_affine, &split)))
            return status;

        *array_size = sizeof(blastn_values_1_1)/sizeof(array_of_8);
        *gap_open_max = 4;
        *gap_extend_max = 2;
    } else if (reward == 3 && penalty == -2) {
        if ((status=s_SplitArrayOf8(blastn_values_3_2, &kValues, &kValues_non_affine, &split)))
            return status;

        *array_size = sizeof(blastn_values_3_2)/sizeof(array_of_8);
        *gap_open_max = 5;
        *gap_extend_max = 5;
    } else if (reward == 4 && penalty == -5) {
        if ((status=s_SplitArrayOf8(blastn_values_4_5, &kValues, &kValues_non_affine, &split)))
            return status;

        *array_size = sizeof(blastn_values_4_5)/sizeof(array_of_8);
        *gap_open_max = 12;
        *gap_extend_max = 8;
    } else if (reward == 5 && penalty == -4) {
        if ((status=s_SplitArrayOf8(blastn_values_5_4, &kValues, &kValues_non_affine, &split)))
            return status;

        *array_size = sizeof(blastn_values_5_4)/sizeof(array_of_8);
        *gap_open_max = 25;
        *gap_extend_max = 10;
    } else  { /* Unsupported reward-penalty */
        status = -1;
        if (error_return) {
            char buffer[256];
            sprintf(buffer, "Substitution scores %d and %d are not supported",
                    reward, penalty);
            Blast_MessageWrite(error_return, eBlastSevError, kBlastMessageNoContext, buffer);
        }
    }
    if (split)
        (*array_size)--;

    if (status == 0)
    {
        if (*array_size > 0)
            *normal = static_cast<array_of_8 *>(BlastMemDup(kValues, (*array_size) * sizeof(array_of_8)));
        if (kValues_non_affine)
            *non_affine = static_cast<array_of_8 *>(BlastMemDup(kValues_non_affine, sizeof(array_of_8)));

        status = s_AdjustGapParametersByGcd(*normal, *non_affine, *array_size, gap_open_max, gap_extend_max, divisor);
    }

    return status;
}

static double
BlastKarlinLtoH(Blast_ScoreFreq* sfp, double lambda)
{
    Int4  score;
    double   H, etonlam, sum, scale;

    double *probs = sfp->sprob;
    Int4 low   = sfp->obs_min,  high  = sfp->obs_max;

    if (lambda < 0.) {
        return -1.;
    }
    if (BlastScoreChk(low, high) != 0) return -1.;

    etonlam = exp( - lambda );
    sum = low * probs[low];
    for( score = low + 1; score <= high; score++ ) {
        sum = score * probs[score] + etonlam * sum;
    }

    scale = BLAST_Powi( etonlam, high );
    if( scale > 0.0 ) {
        H = lambda * sum/scale;
    } else { /* Underflow of exp( -lambda * high ) */
        H = lambda * exp( lambda * high + log(sum) );
    }
    return H;
}

Int2
Blast_KarlinBlkNuclGappedCalc(Blast_KarlinBlk* kbp, Int4 gap_open,
                              Int4 gap_extend, Int4 reward, Int4 penalty,
                              Blast_KarlinBlk* kbp_ungap,
                              bool* round_down,
                              Blast_Message** error_return)
{
    const int kGapOpenIndex = 0;
    const int kGapExtIndex = 1;
    const int kLambdaIndex = 2;
    const int kKIndex = 3;
    const int kHIndex = 4;
    int num_combinations = 0;
    int gap_open_max, gap_extend_max;
    array_of_8* normal=NULL;
    array_of_8* linear=NULL;
    Int2 status = s_GetNuclValuesArray(reward,
                                       penalty,
                                       &num_combinations,
                                       &normal,
                                       &linear,
                                       &gap_open_max,
                                       &gap_extend_max,
                                       round_down,
                                       error_return);

    if (status)
    {
        sfree(normal);
        sfree(linear);
        return status;
    }

            ASSERT(kbp && kbp_ungap);


    /* Try to find the table entry corresponding to input gap costs. */
    if (gap_open == 0 && gap_extend == 0 && linear)
    {
        kbp->Lambda = linear[0][kLambdaIndex];
        kbp->K = linear[0][kKIndex];
        kbp->logK = log(kbp->K);
        kbp->H = linear[0][kHIndex];
    }
    else
    {
        int index=0;
        for (index = 0; index < num_combinations; ++index) {
            if (normal[index][kGapOpenIndex] == gap_open &&
                normal[index][kGapExtIndex] == gap_extend) {
                kbp->Lambda = normal[index][kLambdaIndex];
                kbp->K = normal[index][kKIndex];
                kbp->logK = log(kbp->K);
                kbp->H = normal[index][kHIndex];
                break;
            }
        }

        /* If gap costs are not found in the table, check if they belong to the
        infinite domain, where ungapped values of the parameters can be used. */
        if (index == num_combinations) {
            /* If gap costs are larger than maximal provided in tables, copy
               the values from the ungapped Karlin block. */
            if (gap_open >= gap_open_max && gap_extend >= gap_extend_max) {
                Blast_KarlinBlkCopy(kbp, kbp_ungap);
            } else if (error_return) {
                char buffer[8192];
                int i=0;
                int len=0;
                /* Unsupported gap costs combination. */
                sprintf(buffer, "Gap existence and extension values %ld and %ld "
                                "are not supported for substitution scores %ld and %ld\n",
                        (long) gap_open, (long) gap_extend, (long) reward, (long) penalty);
                for (i = 0; i < num_combinations; ++i)
                {
                    len = (int)strlen(buffer);
                    sprintf(buffer+len, "%ld and %ld are supported existence and extension values\n",
                            (long) normal[i][kGapOpenIndex],  (long) normal[i][kGapExtIndex]);
                }
                len = (int)strlen(buffer);
                sprintf(buffer+len, "%ld and %ld are supported existence and extension values\n",
                        (long) gap_open_max, (long) gap_extend_max);
                len = (int)strlen(buffer);
                sprintf(buffer+len, "Any values more stringent than %ld and %ld are supported\n",
                        (long) gap_open_max, (long) gap_extend_max);
                Blast_MessageWrite(error_return, eBlastSevError, kBlastMessageNoContext, buffer);
                sfree(normal);
                sfree(linear);
                return 1;
            }
        }
    }

    sfree(normal);
    sfree(linear);
    return 0;
}
char* BLAST_StrToUpper(const char* string)
{
    char* retval = NULL;        /* the return value */
    char* p = NULL;             /* auxiliary pointer */

    if ( ! string ) {
        return NULL;
    }

    retval = strdup(string);
    if ( !retval ) {
        return NULL;
    }

    for (p = retval; *p != NULLB; p++) {
        *p = toupper((unsigned char)(*p));
    }
    return retval;
}
#define SAFE_CAST_INT_TO_BOOLEAN(p) (((p) != 0) ? TRUE : FALSE)
bool Blast_SubjectIsNucleotide(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & NUCLEOTIDE_SUBJECT_MASK); }
bool Blast_QueryIsNucleotide(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & NUCLEOTIDE_QUERY_MASK); }

bool Blast_ProgramIsNucleotide(EBlastProgramType p)
{ return Blast_QueryIsNucleotide(p) && Blast_SubjectIsNucleotide(p) ;}


Int2
Blast_ScoreBlkMatrixInit(EBlastProgramType program_number,
                         const BlastScoringOptions* scoring_options,
                         BlastScoreBlk* sbp,
                         GET_MATRIX_PATH get_path)
{
    Int2 status = 0;

    if ( !sbp || !scoring_options ) {

        return 1;
    }

    /* Matrix only scoring is used to disable the greedy extension
       optimisations which avoid use of a full-matrix.  This is
       currently only turned on in RMBlastN -RMH-  */
    sbp->matrix_only_scoring = FALSE;

    if (program_number == eBlastTypeBlastn) {

        BLAST_ScoreSetAmbigRes(sbp, 'N');
        BLAST_ScoreSetAmbigRes(sbp, '-');

        /* If reward/penalty are both zero the calling program is
         * indicating that a matrix must be used to score both the
         * ungapped and gapped alignments.  Set the new
         * matrix_only_scoring.  For now reset reward/penalty to
         * allowed blastn values so that extraneous KA stats can be
         * performed without error. -RMH-
         */
        if ( scoring_options->penalty == 0 && scoring_options->reward == 0 )
        {
            sbp->matrix_only_scoring = TRUE;
            sbp->penalty = BLAST_PENALTY;
            sbp->reward = BLAST_REWARD;
        }else {
            sbp->penalty = scoring_options->penalty;
            sbp->reward = scoring_options->reward;
        }

        if (scoring_options->matrix && *scoring_options->matrix != NULLB) {

            sbp->read_in_matrix = TRUE;
            sbp->name = strdup(scoring_options->matrix);

        } else {
            char buffer[50];
            sbp->read_in_matrix = FALSE;
            sprintf(buffer, "blastn matrix:%ld %ld",
                    (long) sbp->reward, (long) sbp->penalty);
            sbp->name = strdup(buffer);
        }

    } else {
        sbp->read_in_matrix = TRUE;
        BLAST_ScoreSetAmbigRes(sbp, 'X');
        sbp->name = BLAST_StrToUpper(scoring_options->matrix);
    }
    status = Blast_ScoreBlkMatrixFill(sbp, get_path);
    if (status) {
        return status;
    }

    return status;
}
Int2 BlastScoringOptionsSetMatrix(BlastScoringOptions* opts,
                                  const char* matrix_name)
{
    Uint4 i;

    if (matrix_name) {
        sfree(opts->matrix);
        opts->matrix = strdup(matrix_name);
        /* Make it all upper case */
        for (i=0; i<strlen(opts->matrix); ++i)
            opts->matrix[i] = toupper((unsigned char) opts->matrix[i]);
    }
    return 0;
}

Int2
BLAST_FillScoringOptions(BlastScoringOptions* options,
                         EBlastProgramType program_number, bool greedy_extension, Int4 penalty, Int4 reward,
                         const char *matrix, Int4 gap_open, Int4 gap_extend)
{
    if (!options)
        return BLASTERR_INVALIDPARAM;

    if (/*program_number != eBlastTypeBlastn &&
         program_number != eBlastTypePhiBlastn*/
            !Blast_ProgramIsNucleotide(program_number)) {/* protein-protein options. */
        /* If matrix name is not provided, keep the default "BLOSUM62" value filled in
           BlastScoringOptionsNew, otherwise reset it. */
        if (matrix)
            BlastScoringOptionsSetMatrix(options, matrix);
    } else {	/* nucleotide-nucleotide options. */
        if (penalty)
            options->penalty = penalty;
        if (reward)
            options->reward = reward;

        if (greedy_extension) {
            options->gap_open = BLAST_GAP_OPEN_MEGABLAST;
            options->gap_extend = BLAST_GAP_EXTN_MEGABLAST;
        }	else {
            options->gap_open = BLAST_GAP_OPEN_NUCL;
            options->gap_extend = BLAST_GAP_EXTN_NUCL;
        }
    }
    if (gap_open >= 0)
        options->gap_open = gap_open;
    if (gap_extend >= 0)
        options->gap_extend = gap_extend;

    options->program_number = program_number;

    return 0;
}
Int2
BLAST_ScoreSetAmbigRes(BlastScoreBlk* sbp, char ambiguous_res)

{
    Int2 index;
    Uint1* ambig_buffer;

    if (sbp == NULL)
        return 1;

    if (sbp->ambig_occupy >= sbp->ambig_size)
    {
        sbp->ambig_size += 5;
        ambig_buffer = (Uint1 *) calloc(sbp->ambig_size, sizeof(Uint1));
        for (index=0; index<sbp->ambig_occupy; index++)
        {
            ambig_buffer[index] = sbp->ambiguous_res[index];
        }
        sfree(sbp->ambiguous_res);
        sbp->ambiguous_res = ambig_buffer;
    }

    if (sbp->alphabet_code == BLASTAA_SEQ_CODE)
    {
        sbp->ambiguous_res[sbp->ambig_occupy] =
                AMINOACID_TO_NCBISTDAA[toupper((unsigned char) ambiguous_res)];
    }
    else {
        if (sbp->alphabet_code == BLASTNA_SEQ_CODE)
            sbp->ambiguous_res[sbp->ambig_occupy] =
                    IUPACNA_TO_BLASTNA[toupper((unsigned char) ambiguous_res)];
        else if (sbp->alphabet_code == NCBI4NA_SEQ_CODE)
            sbp->ambiguous_res[sbp->ambig_occupy] =
                    IUPACNA_TO_NCBI4NA[toupper((unsigned char) ambiguous_res)];
    }
    (sbp->ambig_occupy)++;


    return 0;
}



Int2
BlastScoringOptionsNew(EBlastProgramType program_number, BlastScoringOptions* *options)
{
*options = (BlastScoringOptions*) calloc(1, sizeof(BlastScoringOptions));

if (*options == NULL)
return BLASTERR_INVALIDPARAM;

if (/*program_number != eBlastTypeBlastn &&
         program_number != eBlastTypePhiBlastn*/
!Blast_ProgramIsNucleotide(program_number)) {/*protein-protein options.*/
(*options)->shift_pen = INT2_MAX;
(*options)->is_ooframe = FALSE;
(*options)->gap_open = BLAST_GAP_OPEN_PROT;
(*options)->gap_extend = BLAST_GAP_EXTN_PROT;
(*options)->matrix = strdup(BLAST_DEFAULT_MATRIX);
} else {	/* nucleotide-nucleotide options. */
(*options)->penalty = BLAST_PENALTY;
(*options)->reward = BLAST_REWARD;
/* This is correct except when greedy extension is used. In that case
   these values would have to be reset. */
(*options)->gap_open = BLAST_GAP_OPEN_NUCL;
(*options)->gap_extend = BLAST_GAP_EXTN_NUCL;
}
if (program_number != eBlastTypeTblastx) {
(*options)->gapped_calculation = TRUE;
}
(*options)->program_number = program_number;
/* By default cross_match-like complexity adjusted scoring is
   turned off.  RMBlastN is currently the only program to use this. -RMH */
(*options)->complexity_adjusted_scoring = FALSE;

return 0;
}



static double
NlmKarlinLambdaNR(double* probs, Int4 d, Int4 low, Int4 high, double lambda0,
                  double tolx, Int4 itmax, Int4 maxNewton, Int4 * itn )
{
    Int4 k;
    double x0, x, a = 0, b = 1;
    double f = 4;  /* Larger than any possible value of the poly in [0,1] */
    Int4 isNewton = 0; /* we haven't yet taken a Newton step. */

    assert( d > 0 );

    x0 = exp( -lambda0 );
    x = ( 0 < x0 && x0 < 1 ) ? x0 : .5;

    for( k = 0; k < itmax; k++ ) { /* all iteration indices k */
        Int4 i;
        double g, fold = f;
        Int4 wasNewton = isNewton; /* If true, then the previous step was a */
        /* Newton step */
        isNewton  = 0;            /* Assume that this step is not */

        /* Horner's rule for evaluating a polynomial and its derivative */
        g = 0;
        f = probs[low];
        for( i = low + d; i < 0; i += d ) {
            g = x * g + f;
            f = f * x + probs[i];
        }
        g = x * g + f;
        f = f * x + probs[0] - 1;
        for( i = d; i <= high; i += d ) {
            g = x * g + f;
            f = f * x + probs[i];
        }
        /* End Horner's rule */

        if( f > 0 ) {
            a = x; /* move the left endpoint */
        } else if( f < 0 ) {
            b = x; /* move the right endpoint */
        } else { /* f == 0 */
            break; /* x is an exact solution */
        }
        if( b - a < 2 * a * ( 1 - b ) * tolx ) {
            /* The midpoint of the interval converged */
            x = (a + b) / 2; break;
        }

        if( k >= maxNewton ||
            /* If convergence of Newton's method appears to be failing; or */
            ( wasNewton && fabs( f ) > .9 * fabs(fold) ) ||
            /* if the previous iteration was a Newton step but didn't decrease
             * f sufficiently; or */
            g >= 0
            /* if a Newton step will move us away from the desired solution */
                ) { /* then */
            /* bisect */
            x = (a + b)/2;
        } else {
            /* try a Newton step */
            double p = - f/g;
            double y = x + p;
            if( y <= a || y >= b ) { /* The proposed iterate is not in (a,b) */
                x = (a + b)/2;
            } else { /* The proposed iterate is in (a,b). Accept it. */
                isNewton = 1;
                x = y;
                if( fabs( p ) < tolx * x * (1-x) ) break; /* Converged */
            } /* else the proposed iterate is in (a,b) */
        } /* else try a Newton step. */
    } /* end for all iteration indices k */
    *itn = k;
    return -log(x)/d;
}

double
Blast_KarlinLambdaNR(Blast_ScoreFreq* sfp, double initialLambdaGuess)
{
    Int4  low;        /* Lowest score (must be negative)  */
    Int4  high;       /* Highest score (must be positive) */
    Int4     itn;
    Int4  i, d;
    double*  sprob;
    double   returnValue;

    low = sfp->obs_min;
    high = sfp->obs_max;
    if (sfp->score_avg >= 0.) {   /* Expected score must be negative */
        return -1.0;
    }
    if (BlastScoreChk(low, high) != 0) return -1.;

    sprob = sfp->sprob;
    /* Find greatest common divisor of all scores */
    for (i = 1, d = -low; i <= high-low && d > 1; ++i) {
        if (sprob[i+low] != 0.0) {
            d = BLAST_Gcd(d, i);
        }
    }
    returnValue =
            NlmKarlinLambdaNR( sprob, d, low, high,
                               initialLambdaGuess,
                               BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
                               20, 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT, &itn );


    return returnValue;
}

Int2
Blast_KarlinBlkUngappedCalc(Blast_KarlinBlk* kbp, Blast_ScoreFreq* sfp)
{


    if (kbp == NULL || sfp == NULL)
        return 1;

    /* Calculate the parameter Lambda */

    kbp->Lambda = Blast_KarlinLambdaNR(sfp, BLAST_KARLIN_LAMBDA0_DEFAULT);
    if (kbp->Lambda < 0.)
        goto ErrExit;


    /* Calculate H */

    kbp->H = BlastKarlinLtoH(sfp, kbp->Lambda);
    if (kbp->H < 0.)
        goto ErrExit;


    /* Calculate K and log(K) */

    kbp->K = BlastKarlinLHtoK(sfp, kbp->Lambda, kbp->H);
    if (kbp->K < 0.)
        goto ErrExit;
    kbp->logK = log(kbp->K);

    /* Normal return */
    return 0;

    ErrExit:
    kbp->Lambda = kbp->H = kbp->K = -1.;
    kbp->logK = HUGE_VAL;
    return 1;
}



Blast_ScoreFreq*
Blast_ScoreFreqNew(Int4 score_min, Int4 score_max)
{
    Blast_ScoreFreq*  sfp;
    Int4  range;

    if (BlastScoreChk(score_min, score_max) != 0)
        return NULL;

    sfp = (Blast_ScoreFreq*) calloc(1, sizeof(Blast_ScoreFreq));
    if (sfp == NULL)
        return NULL;

    range = score_max - score_min + 1;
    sfp->sprob = (double*) calloc(range, sizeof(double));
    if (sfp->sprob == NULL)
    {
        Blast_ScoreFreqFree(sfp);
        return NULL;
    }

    sfp->sprob0 = sfp->sprob;
    sfp->sprob -= score_min;        /* center around 0 */
    sfp->score_min = score_min;
    sfp->score_max = score_max;
    sfp->obs_min = sfp->obs_max = 0;
    sfp->score_avg = 0.0;
    return sfp;
}


static Int2
Blast_ResFreqNormalize(const BlastScoreBlk* sbp, Blast_ResFreq* rfp, double norm)
{
    Int2  alphabet_stop, index;
    double   sum = 0., p;

    if (norm == 0.)
        return 1;

    alphabet_stop = sbp->alphabet_start + sbp->alphabet_size;
    for (index=sbp->alphabet_start; index<alphabet_stop; index++)
    {
        p = rfp->prob[index];
        if (p < 0.)
            return 1;
        sum += p;
    }
    if (sum <= 0.)
        return 0;

    for (index=sbp->alphabet_start; index<alphabet_stop; index++)
    {
        rfp->prob[index] /= sum;
        rfp->prob[index] *= norm;
    }
    return 0;
}


Int2
Blast_ResFreqStdComp(const BlastScoreBlk* sbp, Blast_ResFreq* rfp)
{
    Uint4 index;

    for (index=0; index<DIM(nt_prob); index++) {
        rfp->prob[index] = nt_prob[index].p;
    }



    Blast_ResFreqNormalize(sbp, rfp, 1.0);

    return 0;
}



const int NCBI4NA_TO_BLASTNA[BLASTNA_SIZE] = {
        15,     /* Gap, 0 */
        0,     /* A,   1 */
        1,     /* C,   2 */
        6,     /* M,   3 */
        2,     /* G,   4 */
        4,     /* R,   5 */
        9,     /* S,   6 */
        13,     /* V,   7 */
        3,     /* T,   8 */
        8,     /* W,   9 */
        5,     /* Y,  10 */
        12,     /* H,  11 */
        7,     /* K,  12 */
        11,     /* D,  13 */
        10,     /* B,  14 */
        14      /* N,  15 */
};

const int BLASTNA_TO_NCBI4NA[BLASTNA_SIZE] = {
        1,     /* A, 0 */
        2,     /* C, 1 */
        4,     /* G, 2 */
        8,     /* T, 3 */
        5,     /* R, 4 */
        10,     /* Y, 5 */
        3,     /* M, 6 */
        12,     /* K, 7 */
        9,     /* W, 8 */
        6,     /* S, 9 */
        14,     /* B, 10 */
        13,     /* D, 11 */
        11,     /* H, 12 */
        7,     /* V, 13 */
        15,     /* N, 14 */
        0      /* Gap, 15 */
};

const char BLASTNA_TO_IUPACNA[BLASTNA_SIZE] = {
        'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K',
        'W', 'S', 'B', 'D', 'H', 'V', 'N', '-'
};

const char NCBI4NA_TO_IUPACNA[BLASTNA_SIZE] = {
        '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
        'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
};

const int IUPACNA_TO_BLASTNA[128]={
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15, 0,10, 1,11,15,15, 2,12,15,15, 7,15, 6,14,15,
        15,15, 4, 9, 3,15,13, 8,15, 5,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15};

const int IUPACNA_TO_NCBI4NA[128]={
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
        0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const int AMINOACID_TO_NCBISTDAA[128] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,25, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9,27,10,11,12,13,26,
        14,15,16,17,18,24,19,20,21,22,23, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
        '-','A','B','C','D','E','F','G','H','I','K','L','M',
        'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
        'O', 'J'};
