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

#include "config.h"

/** Joint probabilities for BLOSUM62 */
static constexpr int COMPO_NUM_TRUE_AA = 20;
static const int kReMatrixAdjustmentPseudocounts = 20;
/** relative entropy of BLOSUM62 */
static const double kFixedReBlosum62 = 0.44;

static const double
BLOSUM62_JOINT_PROBS[COMPO_NUM_TRUE_AA][COMPO_NUM_TRUE_AA] =
{
  {2.1497573378347484e-02, 2.3470224274721213e-03, 1.9493235258876179e-03,
   2.1674844853066858e-03, 1.5903351423026848e-03, 1.9242657898716525e-03,
   2.9879059292799641e-03, 5.8158526388051033e-03, 1.1076584657559144e-03,
   3.1880644746334580e-03, 4.4186245468471547e-03, 3.3466571942021082e-03,
   1.3412107617355408e-03, 1.6360627863999076e-03, 2.1568959784943114e-03,
   6.2524987419815400e-03, 3.7180506975672363e-03, 4.0281679108936688e-04,
   1.2999956675626666e-03, 5.0679056444508912e-03},
  {2.3470224274721213e-03, 1.7757465118386322e-02, 1.9786027128591904e-03,
   1.5865480081162602e-03, 3.9365984789376245e-04, 2.4858611089731411e-03,
   2.6933867548771758e-03, 1.7221140903704937e-03, 1.2407382229440791e-03,
   1.2435878276496955e-03, 2.4193952633248727e-03, 6.2339060289407083e-03,
   8.0309461712520876e-04, 9.3181986323789834e-04, 9.5783034332718700e-04,
   2.2660898636037261e-03, 1.7802796534180537e-03, 2.6571979312581875e-04,
   9.2634607111251918e-04, 1.5810185245264004e-03},
  {1.9493235258876179e-03, 1.9786027128591904e-03, 1.4140291972553610e-02,
   3.7201973506001745e-03, 4.3845466068066216e-04, 1.5304436972610567e-03,
   2.2097156829738759e-03, 2.8591871815612977e-03, 1.4301072616183181e-03,
   9.9437221166923172e-04, 1.3690958423974782e-03, 2.4402105140841090e-03,
   5.2943633069226512e-04, 7.5004227978192801e-04, 8.6016459857770028e-04,
   3.1466019144814608e-03, 2.2360795375444384e-03, 1.6159545671597605e-04,
   7.0048422794024819e-04, 1.2014015528772706e-03},
  {2.1674844853066858e-03, 1.5865480081162602e-03, 3.7201973506001745e-03,
   2.1274574617480089e-02, 3.9909227141697264e-04, 1.6481246723433428e-03,
   4.9158017471929655e-03, 2.5221102126636373e-03, 9.5384849402143984e-04,
   1.2347404942429857e-03, 1.5202051791453383e-03, 2.4453087721980561e-03,
   4.6429229320514104e-04, 7.6023722413111566e-04, 1.2373315413524663e-03,
   2.8035127901697272e-03, 1.8961512776990257e-03, 1.6218020183662784e-04,
   5.9842263937853702e-04, 1.3158365660538270e-03},
  {1.5903351423026848e-03, 3.9365984789376245e-04, 4.3845466068066216e-04,
   3.9909227141697264e-04, 1.1931352277704348e-02, 3.0937204045913537e-04,
   3.8338775043186374e-04, 7.6951976030099293e-04, 2.2976387481074697e-04,
   1.0956590131781735e-03, 1.5682982157153873e-03, 5.0124929379033781e-04,
   3.7717165634097634e-04, 5.1389991547056834e-04, 3.6111795849154795e-04,
   1.0432626586831986e-03, 9.3041313726939057e-04, 1.4474923964368156e-04,
   3.4603772624580643e-04, 1.3606607271146112e-03},
  {1.9242657898716525e-03, 2.4858611089731411e-03, 1.5304436972610567e-03,
   1.6481246723433428e-03, 3.0937204045913537e-04, 7.3292255467189687e-03,
   3.5385780499965817e-03, 1.3683038039160171e-03, 1.0489026828741754e-03,
   8.9102936026571569e-04, 1.6174411456311808e-03, 3.0968229715707327e-03,
   7.3993258722701268e-04, 5.4255147972143906e-04, 8.4668181752066874e-04,
   1.8931125300036275e-03, 1.3796838284921874e-03, 2.2737931366728891e-04,
   6.7584155312457842e-04, 1.1660966117775285e-03},
  {2.9879059292799641e-03, 2.6933867548771758e-03, 2.2097156829738759e-03,
   4.9158017471929655e-03, 3.8338775043186374e-04, 3.5385780499965817e-03,
   1.6133927472163669e-02, 1.9380952488713059e-03, 1.3667885452189439e-03,
   1.2192061706431622e-03, 2.0030316026648431e-03, 4.1322603720305197e-03,
   6.7909745467514783e-04, 8.5179405867513139e-04, 1.4216207127018586e-03,
   2.9539180653600089e-03, 2.0493063257644955e-03, 2.6488552587183780e-04,
   8.7044186256788659e-04, 1.6987763526262680e-03},
  {5.8158526388051033e-03, 1.7221140903704937e-03, 2.8591871815612977e-03,
   2.5221102126636373e-03, 7.6951976030099293e-04, 1.3683038039160171e-03,
   1.9380952488713059e-03, 3.7804346453413303e-02, 9.5813607255887238e-04,
   1.3849118546156933e-03, 2.0864716056392773e-03, 2.5392537741810947e-03,
   7.3281559749652399e-04, 1.1976708695723554e-03, 1.3641171883713547e-03,
   3.8342830901664762e-03, 2.1858459940987062e-03, 4.0740829083805248e-04,
   8.3467413018106177e-04, 1.8218235950233687e-03},
  {1.1076584657559144e-03, 1.2407382229440791e-03, 1.4301072616183181e-03,
   9.5384849402143984e-04, 2.2976387481074697e-04, 1.0489026828741754e-03,
   1.3667885452189439e-03, 9.5813607255887238e-04, 9.2802502369336622e-03,
   5.8089627083019206e-04, 9.8696608463236094e-04, 1.1873625842258938e-03,
   3.8264639620910225e-04, 8.1041076335565583e-04, 4.7770135861914477e-04,
   1.1052034635193162e-03, 7.4371746073077327e-04, 1.5168037757411286e-04,
   1.5213771111755425e-03, 6.4882907765797669e-04},
  {3.1880644746334580e-03, 1.2435878276496955e-03, 9.9437221166923172e-04,
   1.2347404942429857e-03, 1.0956590131781735e-03, 8.9102936026571569e-04,
   1.2192061706431622e-03, 1.3849118546156933e-03, 5.8089627083019206e-04,
   1.8441526588740136e-02, 1.1382470627796603e-02, 1.5655862274689192e-03,
   2.5081290988482057e-03, 3.0458868657559346e-03, 1.0068164685944146e-03,
   1.7225081689171561e-03, 2.6953622613315018e-03, 3.6183761166072852e-04,
   1.3821121844492116e-03, 1.1972663837662637e-02},
  {4.4186245468471547e-03, 2.4193952633248727e-03, 1.3690958423974782e-03,
   1.5202051791453383e-03, 1.5682982157153873e-03, 1.6174411456311808e-03,
   2.0030316026648431e-03, 2.0864716056392773e-03, 9.8696608463236094e-04,
   1.1382470627796603e-02, 3.7141460156350926e-02, 2.4634345023228079e-03,
   4.9293545515183088e-03, 5.4151301166464015e-03, 1.4146090399381900e-03,
   2.4277107072013821e-03, 3.3238031308707055e-03, 7.3206640617832933e-04,
   2.2096734692836624e-03, 9.4786263030457313e-03},
  {3.3466571942021082e-03, 6.2339060289407083e-03, 2.4402105140841090e-03,
   2.4453087721980561e-03, 5.0124929379033781e-04, 3.0968229715707327e-03,
   4.1322603720305197e-03, 2.5392537741810947e-03, 1.1873625842258938e-03,
   1.5655862274689192e-03, 2.4634345023228079e-03, 1.6113385590544604e-02,
   9.0876633395557617e-04, 9.4875149773685364e-04, 1.5773020912564391e-03,
   3.1016069999481111e-03, 2.3467014804084987e-03, 2.7198500003555514e-04,
   9.9908866586876396e-04, 1.9360424083099779e-03},
  {1.3412107617355408e-03, 8.0309461712520876e-04, 5.2943633069226512e-04,
   4.6429229320514104e-04, 3.7717165634097634e-04, 7.3993258722701268e-04,
   6.7909745467514783e-04, 7.3281559749652399e-04, 3.8264639620910225e-04,
   2.5081290988482057e-03, 4.9293545515183088e-03, 9.0876633395557617e-04,
   4.0477309321969848e-03, 1.1901770463553603e-03, 4.0824445213456919e-04,
   8.5603787638552766e-04, 1.0095451907679563e-03, 1.9872537223131380e-04,
   5.7145288352831449e-04, 2.3123361470140736e-03},
  {1.6360627863999076e-03, 9.3181986323789834e-04, 7.5004227978192801e-04,
   7.6023722413111566e-04, 5.1389991547056834e-04, 5.4255147972143906e-04,
   8.5179405867513139e-04, 1.1976708695723554e-03, 8.1041076335565583e-04,
   3.0458868657559346e-03, 5.4151301166464015e-03, 9.4875149773685364e-04,
   1.1901770463553603e-03, 1.8277684015431908e-02, 5.2528021756783813e-04,
   1.1939618185901600e-03, 1.1624184369750680e-03, 8.4917468952377874e-04,
   4.2392005745634370e-03, 2.5763052227920180e-03},
  {2.1568959784943114e-03, 9.5783034332718700e-04, 8.6016459857770028e-04,
   1.2373315413524663e-03, 3.6111795849154795e-04, 8.4668181752066874e-04,
   1.4216207127018586e-03, 1.3641171883713547e-03, 4.7770135861914477e-04,
   1.0068164685944146e-03, 1.4146090399381900e-03, 1.5773020912564391e-03,
   4.0824445213456919e-04, 5.2528021756783813e-04, 1.9066033679132538e-02,
   1.6662567934883051e-03, 1.3511005665728870e-03, 1.4152209821874487e-04,
   4.5224391125285910e-04, 1.2451325046931832e-03},
  {6.2524987419815400e-03, 2.2660898636037261e-03, 3.1466019144814608e-03,
   2.8035127901697272e-03, 1.0432626586831986e-03, 1.8931125300036275e-03,
   2.9539180653600089e-03, 3.8342830901664762e-03, 1.1052034635193162e-03,
   1.7225081689171561e-03, 2.4277107072013821e-03, 3.1016069999481111e-03,
   8.5603787638552766e-04, 1.1939618185901600e-03, 1.6662567934883051e-03,
   1.2585947097159817e-02, 4.7004857686835334e-03, 2.8731729176487776e-04,
   1.0299846310599138e-03, 2.3587292053265561e-03},
  {3.7180506975672363e-03, 1.7802796534180537e-03, 2.2360795375444384e-03,
   1.8961512776990257e-03, 9.3041313726939057e-04, 1.3796838284921874e-03,
   2.0493063257644955e-03, 2.1858459940987062e-03, 7.4371746073077327e-04,
   2.6953622613315018e-03, 3.3238031308707055e-03, 2.3467014804084987e-03,
   1.0095451907679563e-03, 1.1624184369750680e-03, 1.3511005665728870e-03,
   4.7004857686835334e-03, 1.2514818886617953e-02, 2.8575770858467209e-04,
   9.4161039895612720e-04, 3.6402328079338207e-03},
  {4.0281679108936688e-04, 2.6571979312581875e-04, 1.6159545671597605e-04,
   1.6218020183662784e-04, 1.4474923964368156e-04, 2.2737931366728891e-04,
   2.6488552587183780e-04, 4.0740829083805248e-04, 1.5168037757411286e-04,
   3.6183761166072852e-04, 7.3206640617832933e-04, 2.7198500003555514e-04,
   1.9872537223131380e-04, 8.4917468952377874e-04, 1.4152209821874487e-04,
   2.8731729176487776e-04, 2.8575770858467209e-04, 6.4699301717154852e-03,
   8.8744160259272527e-04, 3.5578318710317554e-04},
  {1.2999956675626666e-03, 9.2634607111251918e-04, 7.0048422794024819e-04,
   5.9842263937853702e-04, 3.4603772624580643e-04, 6.7584155312457842e-04,
   8.7044186256788659e-04, 8.3467413018106177e-04, 1.5213771111755425e-03,
   1.3821121844492116e-03, 2.2096734692836624e-03, 9.9908866586876396e-04,
   5.7145288352831449e-04, 4.2392005745634370e-03, 4.5224391125285910e-04,
   1.0299846310599138e-03, 9.4161039895612720e-04, 8.8744160259272527e-04,
   1.0246100213822419e-02, 1.5489827890922993e-03},
  {5.0679056444508912e-03, 1.5810185245264004e-03, 1.2014015528772706e-03,
   1.3158365660538270e-03, 1.3606607271146112e-03, 1.1660966117775285e-03,
   1.6987763526262680e-03, 1.8218235950233687e-03, 6.4882907765797669e-04,
   1.1972663837662637e-02, 9.4786263030457313e-03, 1.9360424083099779e-03,
   2.3123361470140736e-03, 2.5763052227920180e-03, 1.2451325046931832e-03,
   2.3587292053265561e-03, 3.6402328079338207e-03, 3.5578318710317554e-04,
   1.5489827890922993e-03, 1.9631915140537640e-02} };

/**
 * @file optimize_target_freq.c
 * Routines for finding an optimal set of target frequencies for the
 * purpose of generating a compositionally adjusted score matrix.  The
 * function for performing this optimization is named
 * OptimizeTargetFrequencies (see below).
 *
 * The optimal target frequencies x minimize the Kullback-Liebler
 * "distance"
 *
 *       sum_k x[k] * ln(x[k]/q[k])
 *
 * from the set q of target frequencies from a standard matrix.  They
 * also satisfy the constraints
 *
 * sum_{i = 0...alphsize - 1} x[i * alphsize + j] = col_sums[j]
 *      for j = 0...alphsize - 1
 *
 * sum_{j = 0...alphsize - 1} x[i * alphsize + j] = row_sums[j]
 *      for i = 1...alphsize - 1
 *
 * where col_sums and row_sums are sets of background frequencies.
 * Note that the second set of equations, the first index for i is one
 * (the i = 0 case holds automatically if the other constraints are
 * satisfied).
 *
 * The coefficient matrix A of these linear equations is used
 * implicitly in this file.  The routines MultiplyByA,
 * MultiplyByAtranspose and ScaledSymmetricProductA all perform
 * operations with A, without explicitly forming the matrix A.  The
 * indices i and j are related to the column indices k of A by the
 * formula
 *
 *    k = i * alphsize + j.
 *
 * Thus A[j][k] = 1 and A[i + alphsize - 1][k] = 1 if i > 0 and
 * A[ell][k] = 0 otherwise.
 *
 * The target frequencies are also usually subject to a constraint on
 * the relative entropy.
 *
 *  relative_entropy =
 *      sum_{ij} x[i * alphsize + j] *
 *                    log(x[i * alphsize + j]/(row_sums[i] * col_sums[j]))
 *
 * but this constraint is optional.
 *
 * REFERENCES
 *
 * Yi-Kuo Yu, John C Wootton, Stephen F. Altschul (2003) The
 * compositional adjustment of amino-acid substitution matrices. Proc
 * Natl Acad Sci USA. 100, 15688-93.
 *
 * Stephen F. Altschul, John C. Wootton, E. Michael Gertz, Richa
 * Agarwala, Aleksandr Morgulis, Alejandro Schaffer and Yi-Kuo Yu
 * (2005) Protein Database Searches Using Compositionally Adjusted
 * Substitution Matrices.  FEBS Journal, 272,5101-9.
 *
 * @author E. Michael Gertz
 */

#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
//#include <algo/blast/core/ncbi_std.h>
#include "../lib/blast/nlm_linear_algebra.h"
#include "cbs.h"
//#include <algo/blast/composition_adjustment/optimize_target_freq.h>

/** bound on error for Newton's method */
/** iteration limit for Newton's method */

 /**
  * Compute the symmetric product A D A^T, where A is the matrix of
  * linear constraints from the problem of generating optimal target
  * frequencies. The matrix A is constant, and is used implicitly in
  * this routine.
  *
  * The symmetric product
  * \f[
  *     A D A^T = \sum_{k = 0}^{n - 1} d_k \times A[:][k] \times A[:][k]^T,
  * \f]
  * where n = alphsize * alphsize and A[:][k] is column k of A.
  *
  * @param alphsize     the size of the alphabet for this minimization
  *                     problem
  * @param W            the product, a matrix of size 2 * alphsize - 1
  * @param diagonal     a vector that represents the diagonal of D, of
  *                     length alphsize * alphsize
  */
static void
ScaledSymmetricProductA(double** W, const double diagonal[], int alphsize)
{
    int rowW, colW;   /* iteration indices over the rows and columns of W */
    int i, j;         /* iteration indices over characters in the alphabet */
    int m;            /* The number of rows in A; also the size of W */

    m = 2 * alphsize - 1;

    for (rowW = 0; rowW < m; rowW++) {
        for (colW = 0; colW <= rowW; colW++) {
            W[rowW][colW] = 0.0;
        }
    }
    for (i = 0; i < alphsize; i++) {
        for (j = 0; j < alphsize; j++) {
            double dd;     /* an individual diagonal element */

            dd = diagonal[i * alphsize + j];

            W[j][j] += dd;
            if (i > 0) {
                W[i + alphsize - 1][j] += dd;
                W[i + alphsize - 1][i + alphsize - 1] += dd;
            }
        }
    }
}


/**
 * Compute the product y = beta * y + alpha * A * x, where A is the
 * constraint matrix for the problem of generating optimal target
 * frequencies.  If beta == 0.0, then y need not be initialized before
 * calling routine.
 *
 * @param beta         a scalar
 * @param alphsize     the alphabet size for this minimization problem
 * @param y            a vector of size 2 * alphsize - 1
 * @param alpha        a scalar
 * @param x            a vector of size alphsize * alphsize
 */
static void
MultiplyByA(double beta, double y[], int alphsize,
    double alpha, const double x[])
{
    int i, j;     /* iteration indices over characters in the alphabet */
    if (beta == 0.0) {
        /* Initialize y to zero, without reading any elements from y */
        for (i = 0; i < 2 * alphsize - 1; i++) {
            y[i] = 0.0;
        }
    }
    else if (beta != 1.0) {
        /* rescale y */
        for (i = 0; i < 2 * alphsize - 1; i++) {
            y[i] *= beta;
        }
    }
    for (i = 0; i < alphsize; i++) {
        for (j = 0; j < alphsize; j++) {
            y[j] += alpha * x[i * alphsize + j];
        }
    }
    for (i = 1; i < alphsize; i++) {
        for (j = 0; j < alphsize; j++) {
            y[i + alphsize - 1] += alpha * x[i * alphsize + j];
        }
    }
}


/**
 * Compute the product y = beta * y + alpha * A^T * x, where A^T is
 * the transpose of the constraint matrix for the problem of
 * generating optimal target frequencies.  If beta == 0.0, then y need
 * not be initialized before calling routine.
 *
 * @param beta         a scalar
 * @param alphsize     the alphabet size for this minimization problem
 * @param y            a vector of size alphsize * alphsize
 * @param alpha        a scalar
 * @param x            a vector of size 2 * alphsize - 1
 */
static void
MultiplyByAtranspose(double beta, double y[], int alphsize,
    double alpha, const double x[])
{
    int i, j;     /* iteration indices over characters in the alphabet */
    int k;        /* index of a row of A transpose (a column of A); also
                      an index into y */

    if (beta == 0.0) {
        /* Initialize y to zero, without reading any elements from y */
        for (k = 0; k < alphsize * alphsize; k++) {
            y[k] = 0.0;
        }
    }
    else if (beta != 1.0) {
        /* rescale y */
        for (k = 0; k < alphsize * alphsize; k++) {
            y[k] *= beta;
        }
    }
    for (i = 0; i < alphsize; i++) {
        for (j = 0; j < alphsize; j++) {
            k = i * alphsize + j;

            y[k] += alpha * x[j];
            if (i > 0) {
                y[k] += alpha * x[i + alphsize - 1];
            }
        }
    }
}


/**
 * Calculate the residuals of the linear constraints of the score
 * matrix optimization problem.
 *
 * @param rA           the residual vector of the linear constraints
 * @param alphsize     the alphabet size for this minimization problem
 * @param x            the substitution probabilities
 * @param row_sums     row sums of the substitution probabilities
 * @param col_sums     column sums of the substitution probabilities
 */
static void
ResidualsLinearConstraints(double rA[], int alphsize, const double x[],
    const double row_sums[], const double col_sums[])
{
    int i;             /* iteration index */

    for (i = 0; i < alphsize; i++) {
        rA[i] = col_sums[i];
    }
    for (i = 1; i < alphsize; i++) {
        rA[i + alphsize - 1] = row_sums[i];
    }
    MultiplyByA(1.0, rA, alphsize, -1.0, x);
}


/**
 * Compute the dual residual vector of the optimization problem.  The
 * dual residual vector is also known as the gradient of the
 * Lagrangian function.
 *
 * @param resids_x     the dual residual vector
 * @param z            dual variables (Lagrange multipliers)
 * @param alphsize     the alphabet size for this optimization problem
 * @param grads        the gradient of the objective function, an
 *                     possibly the nonlinear relative entropy constraint.
 * @param constrain_rel_entropy    if true, then the relative entropy
 *                                 constraint is used in this optimization
 *                                 problem.
 */
static void
DualResiduals(double resids_x[], int alphsize, double** grads,
    const double z[], int constrain_rel_entropy)
{
    int i;                        /* iteration index */
    int n = alphsize * alphsize;  /* size of resids_x */

    if (constrain_rel_entropy) {
        double eta;                /* dual variable for the relative
                                      entropy constraint */
        eta = z[2 * alphsize - 1];
        for (i = 0; i < n; i++) {
            resids_x[i] = -grads[0][i] + eta * grads[1][i];
        }
    }
    else {
        for (i = 0; i < n; i++) {
            resids_x[i] = -grads[0][i];
        }
    }
    MultiplyByAtranspose(1.0, resids_x, alphsize, 1.0, z);
}


/**
 * Calculate the primal and dual residuals the optimization problem,
 * and their combined Euclidean norm.
 *
 * @param rnorm        the Euclidean norm of the residuals
 * @param resids_x     the dual residual vector (gradient of the
 *                     Lagrangian)
 * @param alphsize     size of the alphabet for this minimization
 *                     problem
 * @param resids_z     the primal residual vector (residuals of the
 *                     original constraints)
 * @param values       a vector containing the values of the nonlinear
 *                     functions
 * @param grads        a matrix whose rows are the gradients of the
 *                     nonlinear functions
 * @param row_sums     row sums of the substitution probabilities
 * @param col_sums     column sums of the substitution probabilities
 * @param relative_entropy          target relative entropy
 * @param x            the substitution probabilities
 * @param z            the dual variables (Lagrange multipliers)
 * @param constrain_rel_entropy    if true, the relative entropy
 *                                 constraint is used
 *
 */
static void
CalculateResiduals(double* rnorm,
    double resids_x[],
    int alphsize,
    double resids_z[],
    const double values[],
    double** grads,
    const double row_sums[],
    const double col_sums[],
    const double x[],
    const double z[],
    int constrain_rel_entropy,
    double relative_entropy)
{
    /* Euclidean norms of the primal and dual residuals */
    double norm_resids_z, norm_resids_x;

    DualResiduals(resids_x, alphsize, grads, z, constrain_rel_entropy);
    norm_resids_x = Nlm_EuclideanNorm(resids_x, alphsize * alphsize);

    ResidualsLinearConstraints(resids_z, alphsize, x, row_sums, col_sums);

    if (constrain_rel_entropy) {
        resids_z[2 * alphsize - 1] = relative_entropy - values[1];

        norm_resids_z = Nlm_EuclideanNorm(resids_z, 2 * alphsize);
    }
    else {
        norm_resids_z = Nlm_EuclideanNorm(resids_z, 2 * alphsize - 1);
    }
    *rnorm =
        sqrt(norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z);
}


/**
 * Represents the inverse of the matrix of the linear system that must
 * be solved at every iteration of Newton's method.  If J is the
 * Jacobian of all constraints, then the system has the form
 *
 *     (D     J^T)
 *     (J     0  ),
 *
 * where D is a diagonal matrix.  This matrix may be block reduced to
 * the form.
 *
 *     (D    J^T          )
 *     (0    -J D^{-1} J^T)
 *
 * Neither the matrix nor its inverse are stored explicitly. Rather, a
 * factorization of -J D^{-1} J^T, and the sufficient information to
 * backsolve using this factorization are stored.
 */
typedef struct ReNewtonSystem {
    int alphsize;              /**< the size of the alphabet */
    int constrain_rel_entropy; /**< if true, use the relative entropy
                                    constraint for this optimization
                                    problem */
    double** W;               /**< A lower-triangular matrix
                                    representing a factorization of
                                    the (2,2) block, -J D^{-1} J^T, of
                                    the condensed linear system */
    double* Dinv;             /**< The diagonal elements of the
                                    inverse of the necessarily
                                    diagonal (1,1) block of the linear
                                    system */
    double* grad_re;          /**< the gradient of the
                                    relative-entropy constraint, if
                                    this constraint is used. */
} ReNewtonSystem;


/**
 * Free the memory associated with a ReNewtonSystem.
 *
 * @param newton_system      on entry *newton_system points to the
 *                           system to be freed.  On exit, *newton_system
 *                           is set to NULL.
 */
static void
ReNewtonSystemFree(ReNewtonSystem** newton_system)
{
    if (*newton_system != NULL) {
        Nlm_DenseMatrixFree(&(*newton_system)->W);

        free((*newton_system)->Dinv);
        (*newton_system)->Dinv = NULL;

        free((*newton_system)->grad_re);
        (*newton_system)->grad_re = NULL;

        free(*newton_system);
        *newton_system = NULL;
    }
}


/**
 * Create a new uninitialized ReNewtonSystem; the fields are
 * initialized by the FactorReNewtonSystem procedure.
 * ReNewtonSystemNew and FactorReNewtonSystem are called from only the
 * newt procedure.
 *
 * @param alphsize    the size of the alphabet for this optimization
 *                     problem.
 */
static ReNewtonSystem* ReNewtonSystemNew(int alphsize)
{
    ReNewtonSystem* newton_system;  /* the new ReNewtonSystem */

    newton_system = (ReNewtonSystem*)malloc(sizeof(ReNewtonSystem));
    if (newton_system != NULL) {
        newton_system->alphsize = alphsize;
        newton_system->constrain_rel_entropy = 1;
        newton_system->W = NULL;
        newton_system->Dinv = NULL;
        newton_system->grad_re = NULL;

        newton_system->W = Nlm_LtriangMatrixNew(2 * alphsize);
        if (newton_system->W == NULL)
            goto error_return;
        newton_system->Dinv =
            (double*)malloc(alphsize * alphsize * sizeof(double));
        if (newton_system->Dinv == NULL)
            goto error_return;
        newton_system->grad_re =
            (double*)malloc(alphsize * alphsize * sizeof(double));
        if (newton_system->grad_re == NULL)
            goto error_return;
    }
    goto normal_return;
error_return:
    ReNewtonSystemFree(&newton_system);
normal_return:

    return newton_system;
}


/**
 * Factor the linear system to be solved in this iteration of Newton's
 * method.
 *
 * @param newton_system      holds the factorization
 * @param x                  the primal variables, representing the
 *                           target frequencies.
 * @param z                  the dual variables (Lagrange multipliers)
 * @param grads              gradient of the objective function and
 *                           (if used) the relative entropy constraint
 * @param constrain_rel_entropy    if true, then the relative entropy
 *                                 constraint is used for this optimization
 *                                 problem.
 * @param workspace          an allocated workspace array at least
 *                           as large as the number of target frequencies
 */
static void
FactorReNewtonSystem(ReNewtonSystem* newton_system,
    const double x[],
    const double z[],
    double** grads,
    int constrain_rel_entropy,
    double* workspace)
{
    int i;          /* iteration index */
    int n;          /* the length of x */
    int m;          /* the length of z */

    /* Pointers to fields in newton_systems; the names of the local
     * variables match the names of the fields. */
    double** W = newton_system->W;
    int alphsize = newton_system->alphsize;
    double* Dinv = newton_system->Dinv;
    double* grad_re = newton_system->grad_re;

    n = alphsize * alphsize;
    m = constrain_rel_entropy ? 2 * alphsize : 2 * alphsize - 1;

    newton_system->constrain_rel_entropy = constrain_rel_entropy;

    /* The original system has the form
     *
     *     (D     J^T)
     *     (J     0  ).
     *
     * We block reduce the system to
     *
     *     (D    J^T          )
     *     (0    -J D^{-1} J^T).
     *
     * First we find the inverse of the diagonal matrix D. */

    if (constrain_rel_entropy) {
        double eta;             /* dual variable for the relative
                                   entropy constraint */
        eta = z[m - 1];
        for (i = 0; i < n; i++) {
            Dinv[i] = x[i] / (1 - eta);
        }
    }
    else {
        memcpy(Dinv, x, n * sizeof(double));
    }

    /* Then we compute J D^{-1} J^T; First fill in the part that corresponds
     * to the linear constraints */
    ScaledSymmetricProductA(W, Dinv, alphsize);

    if (constrain_rel_entropy) {
        /* Save the gradient of the relative entropy constraint. */
        memcpy(grad_re, grads[1], n * sizeof(double));

        /* Fill in the part of J D^{-1} J^T that corresponds to the relative
         * entropy constraint. */
        W[m - 1][m - 1] = 0.0;
        for (i = 0; i < n; i++) {
            workspace[i] = Dinv[i] * grad_re[i];

            W[m - 1][m - 1] += grad_re[i] * workspace[i];
        }
        MultiplyByA(0.0, &W[m - 1][0], alphsize, 1.0, workspace);
    }
    /* Factor J D^{-1} J^T and save the result in W. */
    Nlm_FactorLtriangPosDef(W, m);
}


/**
 * Solve the linear system for this iteration of Newton's method, using
 * the matrix factored by the FactorReNewtonSystem routine.
 *
 * @param x               on entry, the dual residuals; on exit, the
 *                        step in the primal variables.
 * @param z               on entry, the primal residuals; on exit, the
 *                        step in the dual variables.
 * @param newton_system   the factored matrix for the Newton system.
 * @param workspace       an allocated workspace array at least
 *                        as large as the number of target frequencies
 */
static void
SolveReNewtonSystem(double x[], double z[],
    const ReNewtonSystem* newton_system, double workspace[])
{
    int i;                     /* iteration index */
    int n;                     /* the size of x */
    int mA;                    /* the number of linear constraints */
    int m;                     /* the size of z */

    /* Local variables that represent fields of newton_system */
    double** W = newton_system->W;
    double* Dinv = newton_system->Dinv;
    double* grad_re = newton_system->grad_re;
    int alphsize = newton_system->alphsize;
    int constrain_rel_entropy = newton_system->constrain_rel_entropy;

    n = alphsize * alphsize;
    mA = 2 * alphsize - 1;
    m = constrain_rel_entropy ? mA + 1 : mA;

    /* Apply the same block reduction to the right-hand side as was
     * applied to the matrix:
     *
     *     rzhat = rz - J D^{-1} rx
     */
    for (i = 0; i < n; i++) {
        workspace[i] = x[i] * Dinv[i];
    }
    MultiplyByA(1.0, z, alphsize, -1.0, workspace);

    if (constrain_rel_entropy) {
        for (i = 0; i < n; i++) {
            z[m - 1] -= grad_re[i] * workspace[i];
        }
    }

    /* Solve for step in z, using the inverse of J D^{-1} J^T */
    Nlm_SolveLtriangPosDef(z, m, W);

    /* Backsolve for the step in x, using the newly-computed step in z.
     *
     *     x = D^{-1) (rx + J\T z)
     */
    if (constrain_rel_entropy) {
        for (i = 0; i < n; i++) {
            x[i] += grad_re[i] * z[m - 1];
        }
    }
    MultiplyByAtranspose(1.0, x, alphsize, 1.0, z);

    for (i = 0; i < n; i++) {
        x[i] *= Dinv[i];
    }
}


/**
 * Evaluate the nonlinear functions and derivatives associated with
 * the problem of generating optimal target frequencies.
 *
 * @param values        values[0] is the value of the objective function
 *                      and values[1] is the relative entropy
 * @param grads         grads[0] is the gradient of the objective function
 *                      and grad[1] is the gradient of the relative entropy
 *                      constraint
 * @param alphsize      the alphabet size for this problem
 * @param x             the primal variables
 * @param q             target frequencies of the standard matrix
 * @param scores        scores as computed using the target frequencies
 *                      of the standard matrix, but the composition
 *                      of the sequences being compared (scaled to
 *                      lambda = 1).
 * @param constrain_rel_entropy     if true, the relative entropy constraint
 *                                  is used in this optimization problem
 */
static void
EvaluateReFunctions(double values[], double** grads, int alphsize,
    const double x[], const double q[],
    const double scores[],
    int constrain_rel_entropy)
{
    int k;         /* iteration index over elements of x, q and scores */
    double temp;   /* holds intermediate values in a computation */

    values[0] = 0.0; values[1] = 0.0;
    for (k = 0; k < alphsize * alphsize; k++) {
        temp = log(x[k] / q[k]);

        values[0] += x[k] * temp;
        grads[0][k] = temp + 1;

        if (constrain_rel_entropy) {
            temp += scores[k];

            values[1] += x[k] * temp;
            grads[1][k] = temp + 1;
        }
    }
}


/**
 * Compute a set of scores (scaled to lambda = 1) from a set of target
 * frequencies and a set of background frequencies.  The target
 * frequencies need not be consistent with the background
 * frequencies.
 *
 * @param scores        the resulting vector of scores, interpreted
 *                      as a matrix stored in row major order
 * @param alphsize      the size of the alphabet
 * @param target_freqs  target frequencies, interpreted as a matrix
 *                      stored in row-major order
 * @param row_freqs     background frequencies of one sequence
 * @param col_freqs     background frequencies of the other sequence
 */
static void
ComputeScoresFromProbs(double scores[],
    int alphsize,
    const double target_freqs[],
    const double row_freqs[],
    const double col_freqs[])
{
    int i, j;     /* iteration indices over characters in the alphabet */
    int k;        /* index into scores and target_freqs */

    for (i = 0; i < alphsize; i++) {
        for (j = 0; j < alphsize; j++) {
            k = i * alphsize + j;

            scores[k] = log(target_freqs[k] / (row_freqs[i] * col_freqs[j]));
        }
    }
}


/* Documented in optimized_target_freq.h */
int
Blast_OptimizeTargetFrequencies(double x[],
    int alphsize,
    int* iterations,
    const double q[],
    const double row_sums[],
    const double col_sums[],
    int constrain_rel_entropy,
    double relative_entropy,
    double tol,
    int maxits)
{
    int its;       /* number of iterations that have been performed */
    int n;         /* number of target frequencies; the size of x */
    int mA;        /* number of linear constraints */
    int m;         /* total number of constraints */

    double         values[2];   /* values of the nonlinear functions
                                   at this iterate */
    double** grads = NULL;     /* gradients of the nonlinear
                                   functions at this iterate */

    ReNewtonSystem*
        newton_system = NULL;   /* factored matrix of the linear
                                   system to be solved at this
                                   iteration */
    double* z = NULL;          /* dual variables (Lagrange multipliers) */
    double* resids_x = NULL;   /* dual residuals (gradient of Lagrangian) */
    double* resids_z = NULL;   /* primal (constraint) residuals */
    double rnorm;               /* norm of the residuals for the
                                   current iterate */
    double* old_scores = NULL; /* a scoring matrix, with lambda = 1,
                                   generated from q, row_sums and
                                   col_sums */
    double* workspace = NULL;  /* A vector for intermediate computations */
    int converged;              /* true if Newton's method converged
                                   to a *minimizer* (strong
                                   second-order point) */
    int status;                 /* the return status */
    n = alphsize * alphsize;
    mA = 2 * alphsize - 1;
    m = constrain_rel_entropy ? mA + 1 : mA;

    newton_system = ReNewtonSystemNew(alphsize);
    if (newton_system == NULL) goto error_return;
    resids_x = (double*)malloc(n * sizeof(double));
    if (resids_x == NULL) goto error_return;
    resids_z = (double*)malloc((mA + 1) * sizeof(double));
    if (resids_z == NULL) goto error_return;
    /* z must be initialized to zero */
    z = (double*)calloc(mA + 1, sizeof(double));
    if (z == NULL) goto error_return;
    old_scores = (double*)malloc(n * sizeof(double));
    if (old_scores == NULL) goto error_return;
    workspace = (double*)malloc(n * sizeof(double));
    if (workspace == NULL) goto error_return;
    grads = Nlm_DenseMatrixNew(2, n);
    if (grads == NULL) goto error_return;

    ComputeScoresFromProbs(old_scores, alphsize, q, row_sums, col_sums);

    /* Use q as the initial value for x */
    memcpy(x, q, n * sizeof(double));
    its = 0;        /* Initialize the iteration count. Note that we may
                       converge in zero iterations if the initial x is
                       optimal. */
    while (its <= maxits) {
        /* Compute the residuals */
        EvaluateReFunctions(values, grads, alphsize, x, q, old_scores,
            constrain_rel_entropy);
        CalculateResiduals(&rnorm, resids_x, alphsize, resids_z, values,
            grads, row_sums, col_sums, x, z,
            constrain_rel_entropy, relative_entropy);

        /* and check convergence; the test correctly handles the case
           in which rnorm is NaN (not a number). */
        if (!(rnorm > tol)) {
            /* We converged at the current iterate */
            break;
        }
        else {
            /* we did not converge, so increment the iteration counter
               and start a new iteration */
            if (++its <= maxits) {
                /* We have not exceeded the maximum number of iterations;
                   take a Newton step. */
                double alpha;       /* a positive number used to scale the
                                       Newton step. */

                FactorReNewtonSystem(newton_system, x, z, grads,
                    constrain_rel_entropy, workspace);
                SolveReNewtonSystem(resids_x, resids_z, newton_system,
                    workspace);

                /* Calculate a value of alpha that ensure that x is
                   positive */
                alpha = Nlm_StepBound(x, n, resids_x, 1.0 / .95);
                alpha *= 0.95;

                Nlm_AddVectors(x, n, alpha, resids_x);
                Nlm_AddVectors(z, m, alpha, resids_z);
            }
        }
    }
    converged = 0;
    if (its <= maxits && rnorm <= tol) {
        /* Newton's iteration converged */
        if (!constrain_rel_entropy || z[m - 1] < 1) {
            /* and the final iterate is a minimizer */
            converged = 1;
        }
    }
    status = converged ? 0 : 1;
    *iterations = its;
    goto normal_return;

error_return:
    status = -1;
    *iterations = 0;
normal_return:

    Nlm_DenseMatrixFree(&grads);
    free(workspace);
    free(old_scores);
    free(z);
    free(resids_z);
    free(resids_x);
    ReNewtonSystemFree(&newton_system);

    return status;
}

/* Documented in composition_adjustment.h. */
void
Blast_ApplyPseudocounts(double* probs20,
    int number_of_observations,
    const double* background_probs20)
{
    int i;                 /* loop index */
    double weight;         /* weight assigned to pseudocounts */
    double sum;            /* sum of the observed frequencies */
    /* pseudocounts as a double */
    double dpseudocounts = kReMatrixAdjustmentPseudocounts;
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
Blast_TrueAaToStdTargetFreqs(double** StdFreq, int StdAlphsize,
    const double* freq)
{
    /* Note I'm using a rough convention for this routine that uppercase
     * letters refer to quantities in the standard (larger) alphabet
     * and lowercase letters refer to the true amino acid (smaller)
     * alphabet.
     */
     /* Shorter names for the sizes of the two alphabets */
    const int small_alphsize = COMPO_NUM_TRUE_AA;
    int A, B;          /* characters in the std (big) alphabet */
    int a, b;          /* characters in the small alphabet */
    double sum;        /* sum of values in target_freq; used to normalize */
    sum = 0.0;
    for (a = 0; a < small_alphsize; a++) {
        for (b = 0; b < small_alphsize; b++) {
            sum += freq[a * 20  + b];
        }
    }
    for (A = 0; A < StdAlphsize; A++) {
        /* for all rows */
        //if (alphaConvert[A] < 0) {
        //    /* the row corresponds to a nonstandard reside */
        //    for (B = 0; B < StdAlphsize; B++) {
        //        StdFreq[A][B] = 0.0;
        //    }
        //}
        //else {
            /* the row corresponds to a standard reside */
            a = A;

            for (B = 0; B < StdAlphsize; B++) {
                /* for all columns */
                if (B < 0) {
                    /* the column corresponds to a nonstandard reside */
                    StdFreq[A][B] = 0.0;
                }
                else {
                    /* the column corresponds to a standard reside */
                    b = B;
                    StdFreq[A][B] = freq[a * 20 + b] / sum;
                }
            }
            /* Set values for two-character ambiguities */
            //StdFreq[A][eBchar] = StdFreq[A][eDchar] + StdFreq[A][eNchar];
            //StdFreq[A][eZchar] = StdFreq[A][eEchar] + StdFreq[A][eQchar];
            //if (StdAlphsize > eJchar) {
            //    StdFreq[A][eJchar] = StdFreq[A][eIchar] + StdFreq[A][eLchar];
            //}
        //}
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
Blast_CalcFreqRatios(double** ratios, int alphsize,
    const double row_prob[], const double col_prob[])
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
s_ScoresStdAlphabet(int** Matrix, int Alphsize,
    const double* target_freq,
    const double row_prob[], const double col_prob[],
    double Lambda)
{
    /* Note: I'm using a rough convention for this routine that uppercase
     * letters refer to quantities in the standard (larger) alphabet
     * and lowercase letters refer to the true amino acid (smaller)
     * alphabet.
     */
    int i;
    /* row and column probabilities in the NCBIstdaa alphabet */
    //double RowProb[COMPO_LARGEST_ALPHABET];
    //double ColProb[COMPO_LARGEST_ALPHABET];
    ///* A double precision score matrix */
    double** Scores = Nlm_DenseMatrixNew(Alphsize, Alphsize);
    if (Scores == NULL) {
        return -1;
    }
    //s_UnpackLetterProbs(RowProb, Alphsize, row_prob);
    //s_SetPairAmbigProbsToSum(RowProb, Alphsize);

    //s_UnpackLetterProbs(ColProb, Alphsize, col_prob);
    //s_SetPairAmbigProbsToSum(ColProb, Alphsize);

    Blast_TrueAaToStdTargetFreqs(Scores, Alphsize, target_freq);
    Blast_CalcFreqRatios(Scores, Alphsize, row_prob, col_prob);
    Blast_FreqRatioToScore(Scores, Alphsize, Alphsize, Lambda);
    //s_SetXUOScores(Scores, Alphsize, RowProb, ColProb);

    s_RoundScoreMatrix(Matrix, Alphsize, Alphsize, Scores);
    Nlm_DenseMatrixFree(&Scores);

    //for (i = 0; i < Alphsize; i++) {
    //    Matrix[i][eStopChar] = StartMatrix[i][eStopChar];
    //    Matrix[eStopChar][i] = StartMatrix[eStopChar][i];
    //}
    return 0;
}

/* Documented in composition_adjustment.h. */
int
Blast_CompositionMatrixAdj(int** matrix,
    int alphsize,
    EMatrixAdjustRule matrix_adjust_rule,
    int length1,
    int length2,
    const double* stdaa_row_probs,
    const double* stdaa_col_probs,
    double lambda)
{
    int iteration_count, status;
    double row_probs[COMPO_NUM_TRUE_AA], col_probs[COMPO_NUM_TRUE_AA];
    /* Target RE when optimizing the matrix; zero if the relative
       entropy should not be constrained. */
    double dummy, desired_re = 0.0;
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
        BLOSUM62_bg);
    Blast_ApplyPseudocounts(col_probs, length2,
        BLOSUM62_bg);

    vector<double> mat_final(20 * 20);

    status =
        Blast_OptimizeTargetFrequencies(mat_final.data(),
            COMPO_NUM_TRUE_AA,
            &iteration_count,
            &BLOSUM62_JOINT_PROBS[0][0],
            row_probs, col_probs,
            (desired_re > 0.0),
            desired_re,
            config.cbs_err_tolerance,
            config.cbs_it_limit);

    if (status != 0)            /* Did not compute the target freqs */
        return status;

    return
        s_ScoresStdAlphabet(matrix, alphsize, mat_final.data(),
            row_probs, col_probs,
            lambda);
}


/** 180 degrees in half a circle */
#define HALF_CIRCLE_DEGREES 180
/** some digits of PI */
#define PI 3.1415926543
/** @{ thresholds used to determine which composition mode to use */
#define QUERY_MATCH_DISTANCE_THRESHOLD 0.16
#define LENGTH_RATIO_THRESHOLD 3.0
#define ANGLE_DEGREE_THRESHOLD 70.0
#define HIGH_PAIR_THRESHOLD 0.4
#define LENGTH_LOWER_THRESHOLD 50
/** @} */


/** Return true if length > 50 and the two most frequent letters
 * occur a total of more that 40% of the time. */
static int
s_HighPairFrequencies(const double* letterProbs, int length)
{
    int i; /*index*/
    double max, second; /*two highest letter probabilities*/

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
s_HighPairEitherSeq(const double* P_query, int length1,
    const double* P_match, int length2)
{
    int result1, result2;

    result1 = s_HighPairFrequencies(P_query, length1);
    result2 = s_HighPairFrequencies(P_match, length2);

    return result1 || result2;
}

/* Documented in composition_adjustment.h. */
double
Blast_GetRelativeEntropy(const double A[], const double B[])
{
    int i;                 /* loop index over letters */
    double temp;           /* intermediate term */
    double value = 0.0;    /* square of relative entropy */

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
    const double* P_query,
    const double* P_match)
{
    EMatrixAdjustRule which_rule; /* which relative entropy mode to
                                     return */
    int i;                       /* loop indices */
    double p_query[COMPO_NUM_TRUE_AA];
    double p_match[COMPO_NUM_TRUE_AA]; /*letter probabilities
                                                for query and match*/
    const double* p_matrix;       /* letter probabilities used in
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

    p_matrix = BLOSUM62_bg; // Blast_GetMatrixBackgroundFreq(matrix_name);

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
        if (//(D_m_q > QUERY_MATCH_DISTANCE_THRESHOLD) &&
            //(len_large / len_small > LENGTH_RATIO_THRESHOLD) &&
            (angle > config.cbs_angle)) {
            which_rule = eCompoScaleOldMatrix;
        }
        else {
            which_rule = eUserSpecifiedRelEntropy;
        }
    }
    return which_rule;
}