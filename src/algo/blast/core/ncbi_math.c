/* $Id: ncbi_math.c 351814 2012-02-01 17:14:10Z ucko $
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
 */

/** @file ncbi_math.c
 * Definitions for portable math library
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: ncbi_math.c 351814 2012-02-01 17:14:10Z ucko $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include "ncbi_math.h"

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

/** size of the next series term that indicates convergence
    in the log and polygamma functions */
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131e-16
#endif

double BLAST_Log1p(double x)
{
   Int4   i;
   double   sum, y;

   if (ABS(x) >= 0.2)
      return log(x+1.);

   /* Limit the loop to 500 terms. */
   for (i=0, sum=0., y=x; i<500 ; ) {
      sum += y/++i;
      if (ABS(y) < DBL_EPSILON)
         break;
      y *= x;
      sum -= y/++i;
      if (y < DBL_EPSILON)
         break;
      y *= x;
   }
   return sum;
}

/** evaluate a specified-order derivative of ln(f(x))
 * @param order The order of derivative to evaluate (0...LOGDERIV_ORDER_MAX).
 *              Derivative order 0 just computes ln(f(x))
 * @param u A list of numerical values of f(x) and its derivatives, all at
 *          the same point x, to be used within the computations 
 * @return 'order'-th derivative of ln(f(x)) or HUGE_VAL if 
 *          order is out of range or u[0] is zero
 */
static double 
s_LogDerivative(Int4 order, double* u)
{
  Int4    i;
  double  y[LOGDERIV_ORDER_MAX+1];
  double  value, tmp;

  if (order < 0 || order > LOGDERIV_ORDER_MAX) {
    return HUGE_VAL;
  }

  if (order > 0 && u[0] == 0.) {
    return HUGE_VAL;
  }
  for (i = 1; i <= order; i++)
    y[i] = u[i] / u[0];

  switch (order) {
  case 0:
    if (u[0] > 0.)
      value = log(u[0]);
    else {
      return HUGE_VAL;
    }
    break;
  case 1:
    value = y[1];
    break;
  case 2:
    value = y[2] - y[1] * y[1];
    break;
  case 3:
    value = y[3] - 3. * y[2] * y[1] + 2. * y[1] * y[1] * y[1];
    break;
  case 4:
    value = y[4] - 4. * y[3] * y[1] - 3. * y[2] * y[2]
      + 12. * y[2] * (tmp = y[1] * y[1]);
    value -= 6. * tmp * tmp;
    break;
  default:
    return HUGE_VAL;
  }
  return value;
}

/** auxiliary values for computation of derivative of ln(gamma(x)) */
static double _default_gamma_coef [] = {
         4.694580336184385e+04,
        -1.560605207784446e+05,
         2.065049568014106e+05,
        -1.388934775095388e+05,
         5.031796415085709e+04,
        -9.601592329182778e+03,
         8.785855930895250e+02,
        -3.155153906098611e+01,
         2.908143421162229e-01,
        -2.319827630494973e-04,
         1.251639670050933e-10
         };

/** Compute a specified-order derivative of ln(gamma(x))
 *  evaluated at some point x 
 *  @param x Value at which derivative will be evaluated
 *  @param order Order of derivative (0...POLYGAMMA_ORDER_MAX)
 *  @return 'order'-th derivative of ln(gamma(x)) at specified x.
 *          Accuracy is to 10 digits for x >= 1
 */
static double
s_GeneralLnGamma(double x, Int4 order)
{
   Int4       i;
   double     xx, tx;
   double     y[POLYGAMMA_ORDER_MAX+1];
   double    tmp, value;
   double    *coef;
   const int      xgamma_dim = DIM(_default_gamma_coef);


   xx = x - 1.;  /* normalize from gamma(x + 1) to xx! */
   tx = xx + xgamma_dim;
   for (i = 0; i <= order; ++i) {
      tmp = tx;
      /* sum the least significant terms first */
      coef = &_default_gamma_coef[xgamma_dim];
      if (i == 0) {
         value = *--coef / tmp;
         while (coef > _default_gamma_coef)
            value += *--coef / --tmp;
      }
      else {
         value = *--coef / BLAST_Powi(tmp, i + 1);
         while (coef > _default_gamma_coef)
            value += *--coef / BLAST_Powi(--tmp, i + 1);
         tmp = BLAST_Factorial(i);
         value *= (i%2 == 0 ? tmp : -tmp);
      }
      y[i] = value;
   }
   ++y[0];
   value = s_LogDerivative(order, y);
   tmp = tx + 0.5;
   switch (order) {
   case 0:
      value += ((NCBIMATH_LNPI+NCBIMATH_LN2) / 2.)
               + (xx + 0.5) * log(tmp) - tmp;
      break;
   case 1:
      value += log(tmp) - xgamma_dim / tmp;
      break;
   case 2:
      value += (tmp + xgamma_dim) / (tmp * tmp);
      break;
   case 3:
      value -= (1. + 2.*xgamma_dim / tmp) / (tmp * tmp);
      break;
   case 4:
      value += 2. * (1. + 3.*xgamma_dim / tmp) / (tmp * tmp * tmp);
      break;
   default:
      tmp = BLAST_Factorial(order - 2) * BLAST_Powi(tmp, 1 - order)
            * (1. + (order - 1) * xgamma_dim / tmp);
      if (order % 2 == 0)
         value += tmp;
      else
         value -= tmp;
      break;
   }
   return value;
}

/** Compute, to 10-digit accuracy, a specified order
 *  derivative of ln(abs(gamma(x))).
 *  @param x value at which derivative will be evaluated
 *  @param order Order of derivative (0...POLYGAMMA_ORDER_MAX)
 *              order = 0, 1, 2, corresponds to ln(gamma), 
 *              digamma, trigamma, etc. Note that the value here
 *              is one less than that suggested by the "di" and "tri" 
 *              prefixes of digamma, trigamma, etc. In other words, 
 *              it is truly the order of the derivative.
 *   @return Computed derivative value, or HUGE_VAL if order is out of range
 */
static double 
s_PolyGamma(double x, Int4 order)
{
   Int4      k;
   double   value, tmp;
   double   y[POLYGAMMA_ORDER_MAX+1], sx;

   if (order < 0 || order > POLYGAMMA_ORDER_MAX) {
      return HUGE_VAL;
   }

   if (x >= 1.)
      return s_GeneralLnGamma(x, order);

   if (x < 0.) {
      value = s_GeneralLnGamma(1. - x, order);
      value = ((order - 1) % 2 == 0 ? value : -value);
      if (order == 0) {
         sx = sin(NCBIMATH_PI * x);
         sx = ABS(sx);
         if ( (x < -0.1 && (ceil(x) == x || sx < 2.*DBL_EPSILON))
               || sx == 0.) {
            return HUGE_VAL;
         }
         value += NCBIMATH_LNPI - log(sx);
      }
      else {
         y[0] = sin(x *= NCBIMATH_PI);
         tmp = 1.;
         for (k = 1; k <= order; k++) {
            tmp *= NCBIMATH_PI;
            y[k] = tmp * sin(x += (NCBIMATH_PI/2.));
         }
         value -= s_LogDerivative(order, y);
      }
   }
   else {
      value = s_GeneralLnGamma(1. + x, order);
      if (order == 0) {
         if (x == 0.) {
            return HUGE_VAL;
         }
         value -= log(x);
      }
      else {
         tmp = BLAST_Factorial(order - 1) * BLAST_Powi(x,  -order);
         value += (order % 2 == 0 ? tmp : - tmp);
      }
   }
   return value;
}

/** Compute ln(abs(gamma(x))) to 10-digit accuracy
 *  @param x Point to evaluate ln(abs(gamma(x)))
 *  @return The function value
 */
static double 
s_LnGamma(double x)
{
   return s_PolyGamma(x, 0);
}
/** Tabulated values of the first few factorials */
static const double kPrecomputedFactorial[] = {
       1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880., 3628800.,
       39916800., 479001600., 6227020800., 87178291200., 1307674368000.,
       20922789888000., 355687428096000., 6402373705728000., 
       121645100408832000., 2432902008176640000., 51090942171709440000., 
       1124000727777607680000., 25852016738884976640000., 
       620448401733239439360000., 15511210043330985984000000.,
       403291461126605635584000000., 10888869450418352160768000000., 
       304888344611713860501504000000., 8841761993739701954543616000000., 
       265252859812191058636308480000000., 8222838654177922817725562880000000.,
       263130836933693530167218012160000000., 
       8683317618811886495518194401280000000., 
       295232799039604140847618609643520000000.
       };
       
double BLAST_Factorial(Int4 n)
{
   if (n < 0)
      return 0.0; /* Undefined! */

   if (n < DIM(kPrecomputedFactorial))
      return kPrecomputedFactorial[n];

   return exp(s_LnGamma((double)(n + 1)));
}

double BLAST_LnGammaInt(Int4 n)
{
   if ( (n > 1) && (n < DIM(kPrecomputedFactorial) ) ) {
      return log(kPrecomputedFactorial[n-1]);
   }

   return s_LnGamma((double)n);
}

/*
   Romberg numerical integrator
   Reference:
      Francis Scheid (1968)
      Schaum's Outline Series
      Numerical Analysis, p. 115
      McGraw-Hill Book Company, New York
*/

/** Make a parametrized function appear to have only one variable */
#define F(x)  ((*f)((x), fargs))
/** Maximum number of diagonals in the Romberg array */
#define MAX_DIAGS 20

double BLAST_RombergIntegrate(double (*f) (double,void*), void* fargs, double p, double q, double eps, Int4 epsit, Int4 itmin)

{
   double   romb[MAX_DIAGS];   /* present list of Romberg values */
   double   h;   /* mesh-size */
   Int4      i, j, k, npts;
   long   n;   /* 4^(error order in romb[i]) */
   Int4      epsit_cnt = 0, epsck;
   double   y;
   double   x;
   double   sum;

   /* itmin = min. no. of iterations to perform */
   itmin = MAX(1, itmin);
   itmin = MIN(itmin, MAX_DIAGS-1);

   /* epsit = min. no. of consecutive iterations that must satisfy epsilon */
   epsit = MAX(epsit, 1); /* default = 1 */
   epsit = MIN(epsit, 3); /* if > 3, the problem needs more prior analysis */

   epsck = itmin - epsit;

   npts = 1;
   h = q - p;
   x = F(p);
   if (ABS(x) == HUGE_VAL)
      return x;
   y = F(q);
   if (ABS(y) == HUGE_VAL)
      return y;
   romb[0] = 0.5 * h * (x + y);   /* trapezoidal rule */
   for (i = 1; i < MAX_DIAGS; ++i, npts *= 2, h *= 0.5) {
      sum = 0.;   /* sum of ordinates for 
                     x = p+0.5*h, p+1.5*h, ..., q-0.5*h */
      for (k = 0, x = p+0.5*h; k < npts; k++, x += h) {
         y = F(x);
         if (ABS(y) == HUGE_VAL)
            return y;
         sum += y;
      }
      romb[i] = 0.5 * (romb[i-1] + h*sum); /* new trapezoidal estimate */

      /* Update Romberg array with new column */
      for (n = 4, j = i-1; j >= 0; n *= 4, --j)
         romb[j] = (n*romb[j+1] - romb[j]) / (n-1);

      if (i > epsck) {
         if (ABS(romb[1] - romb[0]) > eps * ABS(romb[0])) {
            epsit_cnt = 0;
            continue;
         }
         ++epsit_cnt;
         if (i >= itmin && epsit_cnt >= epsit)
            return romb[0];
      }
   }
   return HUGE_VAL;
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

Int4 
BLAST_Gdb3(Int4* a, Int4* b, Int4* c)
{
    Int4 g;
    if (*b == 0) 
        g = BLAST_Gcd(*a, *c);
    else 
        g = BLAST_Gcd(*a, BLAST_Gcd(*b, *c));
    if (g > 1) {
        *a /= g;
        *b /= g;
        *c /= g;
    }
    return g;
}

long BLAST_Nint(double x)
{
   x += (x >= 0. ? 0.5 : -0.5);
   return (long)x;
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

double BLAST_LnFactorial (double x) {

    if( x <= 0.0)
        return 0.0;
    else
        return s_LnGamma(x + 1.0);
        
}

#ifdef _NCBILCL_ /* C Toolkit */
#  ifdef OS_UNIX
#    define HAVE_ERF 1
#  endif
#  ifdef IS_BIG_ENDIAN
#    define WORDS_BIGENDIAN 1
#  endif
#endif
/* erf/erfc implementation, cloned from the C++ Toolkit's util directory. */
#define NCBI_Erf  BLAST_Erf
#define NCBI_ErfC BLAST_ErfC
#define NCBI_INCLUDE_NCBI_ERF_C 1
#include "ncbi_erf.c"
