/* $Id: ncbi_math.h 351792 2012-02-01 16:20:19Z ucko $
 * ===========================================================================
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
 * ===========================================================================
 *
 * Authors:  Gish, Kans, Ostell, Schuler
 *
 * Version Creation Date:   10/23/91
 *
 * ==========================================================================
 */

/** @file ncbi_math.h
 * Prototypes for portable math library (ported from C Toolkit)
 */

#ifndef ALGO_BLAST_CORE__NCBIMATH
#define ALGO_BLAST_CORE__NCBIMATH

/** Natural logarithm with shifted input
 *  @param x input operand (x > -1)
 *  @return log(x+1)
 */
double BLAST_Log1p (double x);

/** Exponentional with base e 
 *  @param x input operand
 *  @return exp(x) - 1
 */
double BLAST_Expm1 (double x);

/** Factorial function
 *  @param n input operand
 *  @return (double)(1 * 2 * 3 * ... * n)
 */
double BLAST_Factorial(int32_t n);

/** Logarithm of the factorial 
 *  @param x input operand
 *  @return log(1 * 2 * 3 * ... * x)
 */
double BLAST_LnFactorial (double x);

/** log(gamma(n)), integral n 
 *  @param n input operand
 *  @return log(1 * 2 * 3 * ... (n-1))
 */
double BLAST_LnGammaInt (int32_t n);

/** Romberg numerical integrator 
 *  @param f Pointer to the function to integrate; the first argument
 *               is the variable to integrate over, the second is a pointer
 *               to a list of additional arguments that f may need
 *  @param fargs Pointer to an array of extra arguments or parameters
 *               needed to compute the function to be integrated. None
 *               of the items in this list may vary over the region
 *               of integration
 *  @param p Left-hand endpoint of the integration interval
 *  @param q Right-hand endpoint of the integration interval
 *           (q is assumed > p)
 *  @param eps The relative error tolerance that indicates convergence
 *  @param epsit The number of consecutive diagonal entries in the 
 *               Romberg array whose relative difference must be less than
 *               eps before convergence is assumed. This is presently 
 *               limited to 1, 2, or 3
 *  @param itmin The minimum number of diagnonal Romberg entries that
 *               will be computed
 *  @return The computed integral of f() between p and q
 */
double BLAST_RombergIntegrate (double (*f) (double, void*), 
                               void* fargs, double p, double q, 
                               double eps, int32_t epsit, int32_t itmin);

/** Greatest common divisor 
 *  @param a First operand (any integer)
 *  @param b Second operand (any integer)
 *  @return The largest integer that evenly divides a and b
 */
int32_t BLAST_Gcd (int32_t a, int32_t b);

/** Divide 3 numbers by their greatest common divisor
 * @param a First integer [in] [out]
 * @param b Second integer [in] [out]
 * @param c Third integer [in] [out]
 * @return The greatest common divisor
 */
int32_t BLAST_Gdb3(int32_t* a, int32_t* b, int32_t* c);

/** Nearest integer 
 *  @param x Input to round (rounded value must be representable
 *           as a 32-bit signed integer)
 *  @return floor(x + 0.5);
 */
long BLAST_Nint (double x);

/** Integral power of x 
 * @param x floating-point base of the exponential
 * @param n (integer) exponent
 * @return x multiplied by itself n times
 */
double BLAST_Powi (double x, int32_t n);

/** The error function of x: the integral from 0 to x of e(-t*t) dt,
 *  scaled by 2/sqrt(pi) to fall within the range (-1,1). */
double BLAST_Erf (double x);

/** The complementary error function of x: 1 - erf(x), but calculated
 *  more accurately for large x (where erf(x) approaches unity). */
double BLAST_ErfC (double x);

/** Number of derivatives of log(x) to carry in gamma-related 
    computations */
#define LOGDERIV_ORDER_MAX	4  
/** Number of derivatives of polygamma(x) to carry in gamma-related 
    computations for non-integral values of x */
#define POLYGAMMA_ORDER_MAX	LOGDERIV_ORDER_MAX

/** value of pi is only used in gamma-related computations */
#define NCBIMATH_PI	3.1415926535897932384626433832795

/** Natural log(2) */
#define NCBIMATH_LN2	0.69314718055994530941723212145818
/** Natural log(PI) */
#define NCBIMATH_LNPI	1.1447298858494001741434273513531

#endif /* !ALGO_BLAST_CORE__NCBIMATH */
