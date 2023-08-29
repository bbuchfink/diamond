#pragma once
#define _USE_MATH_DEFINES
#include <math.h>

#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */

enum {
    GSL_SUCCESS = 0,
    GSL_FAILURE = -1,
    GSL_CONTINUE = -2,  /* iteration has not converged */
    GSL_EDOM = 1,   /* input domain error, e.g sqrt(-1) */
    GSL_ERANGE = 2,   /* output range error, e.g. exp(1e100) */
    GSL_EFAULT = 3,   /* invalid pointer */
    GSL_EINVAL = 4,   /* invalid argument supplied by user */
    GSL_EFAILED = 5,   /* generic failure */
    GSL_EFACTOR = 6,   /* factorization failed */
    GSL_ESANITY = 7,   /* sanity check failed - shouldn't happen */
    GSL_ENOMEM = 8,   /* malloc failed */
    GSL_EBADFUNC = 9,   /* problem with user-supplied function */
    GSL_ERUNAWAY = 10,  /* iterative process is out of control */
    GSL_EMAXITER = 11,  /* exceeded max number of iterations */
    GSL_EZERODIV = 12,  /* tried to divide by zero */
    GSL_EBADTOL = 13,  /* user specified an invalid tolerance */
    GSL_ETOL = 14,  /* failed to reach the specified tolerance */
    GSL_EUNDRFLW = 15,  /* underflow */
    GSL_EOVRFLW = 16,  /* overflow  */
    GSL_ELOSS = 17,  /* loss of accuracy */
    GSL_EROUND = 18,  /* failed because of roundoff error */
    GSL_EBADLEN = 19,  /* matrix, vector lengths are not conformant */
    GSL_ENOTSQR = 20,  /* matrix not square */
    GSL_ESING = 21,  /* apparent singularity detected */
    GSL_EDIVERGE = 22,  /* integral or series is divergent */
    GSL_EUNSUP = 23,  /* requested feature is not supported by the hardware */
    GSL_EUNIMPL = 24,  /* requested feature not (yet) implemented */
    GSL_ECACHE = 25,  /* cache limit exceeded */
    GSL_ETABLE = 26,  /* table limit exceeded */
    GSL_ENOPROG = 27,  /* iteration is not making progress towards solution */
    GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
    GSL_ETOLF = 29,  /* cannot reach the specified tolerance in F */
    GSL_ETOLX = 30,  /* cannot reach the specified tolerance in X */
    GSL_ETOLG = 31,  /* cannot reach the specified tolerance in gradient */
    GSL_EOF = 32   /* end of file */
};

void
gsl_stream_printf(const char* label, const char* file, int line,
    const char* reason)
{
    fprintf(stderr, "gsl: %s:%d: %s: %s\n", file, line, label, reason);
}

void
gsl_error(const char* reason, const char* file, int line, int gsl_errno)
{
    gsl_stream_printf("ERROR", file, line, reason);

    fflush(stdout);
    fprintf(stderr, "Default GSL error handler invoked.\n");
    fflush(stderr);

    abort();
}

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

#define EVAL_RESULT(fn) \
   gsl_sf_result result; \
   int status = fn; \
   if (status != GSL_SUCCESS) { \
     GSL_ERROR_VAL(#fn, status, result.val); \
   } ; \
   return result.val;

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

struct gsl_sf_result_struct {
    double val;
    double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
    double val;
    double err;
    int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;

__END_DECLS

struct cheb_series_struct {
    double* c;   /* coefficients                */
    int order;    /* order of expansion          */
    double a;     /* lower interval point        */
    double b;     /* upper interval point        */
    int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;

static inline int
cheb_eval_e(const cheb_series* cs,
    const double x,
    gsl_sf_result* result)
{
    int j;
    double d = 0.0;
    double dd = 0.0;

    double y = (2.0 * x - cs->a - cs->b) / (cs->b - cs->a);
    double y2 = 2.0 * y;

    double e = 0.0;

    for (j = cs->order; j >= 1; j--) {
        double temp = d;
        d = y2 * d - dd + cs->c[j];
        e += fabs(y2 * temp) + fabs(dd) + fabs(cs->c[j]);
        dd = temp;
    }

    {
        double temp = d;
        d = y * d - dd + 0.5 * cs->c[0];
        e += fabs(y * temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
    }

    result->val = d;
    result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

    return GSL_SUCCESS;
}

/* Chebyshev fit for erfc((t+1)/2), -1 < t < 1
 */
static double erfc_xlt1_data[20] = {
  1.06073416421769980345174155056,
 -0.42582445804381043569204735291,
  0.04955262679620434040357683080,
  0.00449293488768382749558001242,
 -0.00129194104658496953494224761,
 -0.00001836389292149396270416979,
  0.00002211114704099526291538556,
 -5.23337485234257134673693179020e-7,
 -2.78184788833537885382530989578e-7,
  1.41158092748813114560316684249e-8,
  2.72571296330561699984539141865e-9,
 -2.06343904872070629406401492476e-10,
 -2.14273991996785367924201401812e-11,
  2.22990255539358204580285098119e-12,
  1.36250074650698280575807934155e-13,
 -1.95144010922293091898995913038e-14,
 -6.85627169231704599442806370690e-16,
  1.44506492869699938239521607493e-16,
  2.45935306460536488037576200030e-18,
 -9.29599561220523396007359328540e-19
};
static cheb_series erfc_xlt1_cs = {
  erfc_xlt1_data,
  19,
  -1, 1,
  12
};

/* Chebyshev fit for erfc(x) exp(x^2), 1 < x < 5, x = 2t + 3, -1 < t < 1
 */
static double erfc_x15_data[25] = {
  0.44045832024338111077637466616,
 -0.143958836762168335790826895326,
  0.044786499817939267247056666937,
 -0.013343124200271211203618353102,
  0.003824682739750469767692372556,
 -0.001058699227195126547306482530,
  0.000283859419210073742736310108,
 -0.000073906170662206760483959432,
  0.000018725312521489179015872934,
 -4.62530981164919445131297264430e-6,
  1.11558657244432857487884006422e-6,
 -2.63098662650834130067808832725e-7,
  6.07462122724551777372119408710e-8,
 -1.37460865539865444777251011793e-8,
  3.05157051905475145520096717210e-9,
 -6.65174789720310713757307724790e-10,
  1.42483346273207784489792999706e-10,
 -3.00141127395323902092018744545e-11,
  6.22171792645348091472914001250e-12,
 -1.26994639225668496876152836555e-12,
  2.55385883033257575402681845385e-13,
 -5.06258237507038698392265499770e-14,
  9.89705409478327321641264227110e-15,
 -1.90685978789192181051961024995e-15,
  3.50826648032737849245113757340e-16
};
static cheb_series erfc_x15_cs = {
  erfc_x15_data,
  24,
  -1, 1,
  16
};

/* Chebyshev fit for erfc(x) x exp(x^2), 5 < x < 10, x = (5t + 15)/2, -1 < t < 1
 */
static double erfc_x510_data[20] = {
  1.11684990123545698684297865808,
  0.003736240359381998520654927536,
 -0.000916623948045470238763619870,
  0.000199094325044940833965078819,
 -0.000040276384918650072591781859,
  7.76515264697061049477127605790e-6,
 -1.44464794206689070402099225301e-6,
  2.61311930343463958393485241947e-7,
 -4.61833026634844152345304095560e-8,
  8.00253111512943601598732144340e-9,
 -1.36291114862793031395712122089e-9,
  2.28570483090160869607683087722e-10,
 -3.78022521563251805044056974560e-11,
  6.17253683874528285729910462130e-12,
 -9.96019290955316888445830597430e-13,
  1.58953143706980770269506726000e-13,
 -2.51045971047162509999527428316e-14,
  3.92607828989125810013581287560e-15,
 -6.07970619384160374392535453420e-16,
  9.12600607264794717315507477670e-17
};
static cheb_series erfc_x510_cs = {
  erfc_x510_data,
  19,
  -1, 1,
  12
};

static double erfc8_sum(double x)
{
    /* estimates erfc(x) valid for 8 < x < 100 */
    /* This is based on index 5725 in Hart et al */

    static double P[] = {
        2.97886562639399288862,
        7.409740605964741794425,
        6.1602098531096305440906,
        5.019049726784267463450058,
        1.275366644729965952479585264,
        0.5641895835477550741253201704
    };
    static double Q[] = {
        3.3690752069827527677,
        9.608965327192787870698,
        17.08144074746600431571095,
        12.0489519278551290360340491,
        9.396034016235054150430579648,
        2.260528520767326969591866945,
        1.0
    };
    double num = 0.0, den = 0.0;
    int i;

    num = P[5];
    for (i = 4; i >= 0; --i) {
        num = x * num + P[i];
    }
    den = Q[6];
    for (i = 5; i >= 0; --i) {
        den = x * den + Q[i];
    }

    return num / den;
}

inline
static double erfc8(double x)
{
    double e;
    e = erfc8_sum(x);
    e *= exp(-x * x);
    return e;
}

inline
static double log_erfc8(double x)
{
    double e;
    e = erfc8_sum(x);
    e = log(e) - x * x;
    return e;
}

int gsl_sf_erfc_e(double x, gsl_sf_result* result)
{
    const double ax = fabs(x);
    double e_val, e_err;

    /* CHECK_POINTER(result) */

    if (ax <= 1.0) {
        double t = 2.0 * ax - 1.0;
        gsl_sf_result c;
        cheb_eval_e(&erfc_xlt1_cs, t, &c);
        e_val = c.val;
        e_err = c.err;
    }
    else if (ax <= 5.0) {
        double ex2 = exp(-x * x);
        double t = 0.5 * (ax - 3.0);
        gsl_sf_result c;
        cheb_eval_e(&erfc_x15_cs, t, &c);
        e_val = ex2 * c.val;
        e_err = ex2 * (c.err + 2.0 * fabs(x) * GSL_DBL_EPSILON);
    }
    else if (ax < 10.0) {
        double exterm = exp(-x * x) / ax;
        double t = (2.0 * ax - 15.0) / 5.0;
        gsl_sf_result c;
        cheb_eval_e(&erfc_x510_cs, t, &c);
        e_val = exterm * c.val;
        e_err = exterm * (c.err + 2.0 * fabs(x) * GSL_DBL_EPSILON + GSL_DBL_EPSILON);
    }
    else {
        e_val = erfc8(ax);
        e_err = (x * x + 1.0) * GSL_DBL_EPSILON * fabs(e_val);
    }

    if (x < 0.0) {
        result->val = 2.0 - e_val;
        result->err = e_err;
        result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
    else {
        result->val = e_val;
        result->err = e_err;
        result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }

    return GSL_SUCCESS;
}

inline int gsl_sf_log_erfc_e(double x, gsl_sf_result* result)
{
    /* CHECK_POINTER(result) */

    if (x * x < 10.0 * GSL_ROOT6_DBL_EPSILON) {
        const double y = x / M_SQRTPI;
        /* series for -1/2 Log[Erfc[Sqrt[Pi] y]] */
        const double c3 = (4.0 - M_PI) / 3.0;
        const double c4 = 2.0 * (1.0 - M_PI / 3.0);
        const double c5 = -0.001829764677455021;  /* (96.0 - 40.0*M_PI + 3.0*M_PI*M_PI)/30.0  */
        const double c6 = 0.02629651521057465;   /* 2.0*(120.0 - 60.0*M_PI + 7.0*M_PI*M_PI)/45.0 */
        const double c7 = -0.01621575378835404;
        const double c8 = 0.00125993961762116;
        const double c9 = 0.00556964649138;
        const double c10 = -0.0045563339802;
        const double c11 = 0.0009461589032;
        const double c12 = 0.0013200243174;
        const double c13 = -0.00142906;
        const double c14 = 0.00048204;
        double series = c8 + y * (c9 + y * (c10 + y * (c11 + y * (c12 + y * (c13 + c14 * y)))));
        series = y * (1.0 + y * (1.0 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * (c7 + y * series)))))));
        result->val = -2.0 * series;
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        return GSL_SUCCESS;
    }
    /*
    don't like use of log1p(); added above series stuff for small x instead, should be ok [GJ]
    else if (fabs(x) < 1.0) {
      gsl_sf_result result_erf;
      gsl_sf_erf_e(x, &result_erf);
      result->val  = log1p(-result_erf.val);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    */
    else if (x > 8.0) {
        result->val = log_erfc8(x);
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        return GSL_SUCCESS;
    }
    else {
        gsl_sf_result result_erfc;
        gsl_sf_erfc_e(x, &result_erfc);
        result->val = log(result_erfc.val);
        result->err = fabs(result_erfc.err / result_erfc.val);
        result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        return GSL_SUCCESS;
    }
}

inline double gsl_sf_log_erfc(double x)
{
    EVAL_RESULT(gsl_sf_log_erfc_e(x, &result));
}