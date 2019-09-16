// Copyright 2015 Yutaro Konta

#include "LambdaCalculator.hh"
#include <vector>
#include <cassert>
#include <cstdio>  // sprintf
#include <cstdlib>
#include <cfloat>
#include <cmath>
//#include <iostream>
using namespace std;

static double roundToFewDigits(double x)
{
  // This rounding fixes some inaccuracies, e.g. for DNA with a simple
  // match/mismatch matrix it's likely to make all the probabilities
  // exactly 0.25, as they should be.
  char buffer[32];
  sprintf(buffer, "%g", x);
  return atof(buffer);
}

static double** makematrix(int m, int n, double val)
{
  double** mat = new double* [m];
  for (int i=0; i<m; i++)
    {
      mat[i] = new double [n];
      for (int j=0; j<n; j++)
        mat[i][j] = val;
    }
  return mat;
}

static void deletematrix(double** a, int m)
{
  for (int i=0; i<m; i++)
    delete[] a[i];
  delete[] a;
}

static double summatrix(double** a, int m, int n)
{
  double s = 0;
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      s += a[i][j];
  return s;
}

static int max_index(double** a, int n, int i)
{
  double max = -DBL_MAX;
  int maxindex = -1;

  for (int j=i; j<n; j++)
    {
      if (fabs(a[j][i]) > max)
        {
          max = fabs(a[j][i]);
          maxindex = j;
        }
    }
  return maxindex;
}

static void swap_matrix(double** a, int i, int j)
{
  double* v = a[i];
  a[i] = a[j];
  a[j] = v;
}

static void swap_vector(int* a, int i, int j)
{
  int v = a[i];
  a[i] = a[j];
  a[j] = v;
}

static bool lu_pivoting(double** a, int* idx, int n)
{
  int v;

  for (int i=0; i<n; i++)
    idx[i] = i;

  for (int i=0; i<n; i++)
    {
      v = max_index(a, n, i);
      assert(v >= 0);
      if (fabs(a[v][i]) < 1e-10)
        {
          return false; // singular matrix!
        }

      swap_matrix(a, i, v);
      swap_vector(idx, i, v);

      a[i][i] = 1.0/a[i][i];
      for (int j=i+1; j<n; j++)
        {
          a[j][i] = a[j][i] * a[i][i];
          for (int k=i+1; k<n; k++)
            a[j][k] = a[j][k] - a[j][i] * a[i][k];
        }
    }
  return true;
}

static void solvep(double **a, double *x, double *b, int n)
{
  double *y = new double [n];

  for (int i=0; i<n; i++)
    {
      y[i] = b[i];
      for (int j=0; j<i; j++)
        y[i] -= a[i][j] * y[j];
    }

  for (int i=n-1; i>=0; i--)
    {
      x[i] = y[i];
      for (int j=i+1; j<n; j++)
        x[i] -= a[i][j] * x[j];
      x[i] *= a[i][i]; // needed because a[i][i] is inverted
    }
  delete[] y;
}

static void transpose(double **a, int n)
{
  double v;
  for (int i=0; i<n; i++)
    {
      for (int j=0; j<i; j++)
        {
          v = a[i][j];
          a[i][j] = a[j][i];
          a[j][i] = v;
        }
    }
}

static bool invert(double **a, double **inv, int n)
{
  int* idx = new int [n];

  double **e = makematrix(n,n,0);

  if(!lu_pivoting(a, idx, n))
    return false;

  for (int i=0; i<n; i++)
    e[idx[i]][i] = 1; // transposed

  delete[] idx;

  for (int i=0; i<n; i++)
    solvep(a, inv[i], e[i], n);

  deletematrix(e, n);
  transpose(inv, n); // transpose inv to make the inverted matrix of a
  return true;
}

static bool calculate_inv_sum(double **matrix, int alpha_size, double tau, double* inv_sum)
{
  double **m = makematrix(alpha_size, alpha_size, 0);
  double **y = makematrix(alpha_size, alpha_size, 0);

  for (int i=0; i<alpha_size; i++)
    for (int j=0; j<alpha_size; j++)
      m[i][j] = exp(tau * matrix[i][j]);

  if(!invert(m, y, alpha_size))
    return false;

  *inv_sum = summatrix(y, alpha_size, alpha_size);

  deletematrix(m, alpha_size);
  deletematrix(y, alpha_size);
  return true;
}

namespace cbrc{

void LambdaCalculator::setBad(){
  lambda_ = -1;
  letterProbs1_.clear();
  letterProbs2_.clear();
}

bool LambdaCalculator::find_ub(double **matrix, int alpha_size, double *ub)
{
  double r_max_min = DBL_MAX;
  double c_max_min = DBL_MAX;
  double r_max;
  double c_max;

  double r_min;
  double c_min;

  int l_r = 0;
  int l_c = 0;

  for (int i=0; i<alpha_size; i++)
    {
      r_max = -DBL_MAX;
      r_min = DBL_MAX;
      for (int j=0; j<alpha_size; j++)
        {
          if (matrix[i][j] > r_max)
            r_max = matrix[i][j];

          if (matrix[i][j] < r_min)
            r_min = matrix[i][j];
        }
      if (r_max == 0 && r_min == 0)
        l_r++;
      else if (r_max <= 0 || r_min >= 0)
        return false;
      else if (r_max < r_max_min)
        r_max_min = r_max;
    }

  for (int j=0; j<alpha_size; j++)
    {
      c_max = -DBL_MAX;
      c_min = DBL_MAX;
      for (int i=0; i<alpha_size; i++)
        {
          if (matrix[i][j] > c_max)
            c_max = matrix[i][j];

          if (matrix[i][j] < c_min)
            c_min = matrix[i][j];
        }
      if (c_max == 0 && c_min == 0)
        l_c++;
      else if (c_max <= 0 || c_min >= 0)
        return false;
      else if (c_max < c_max_min)
        c_max_min = c_max;
    }

  if (l_r == alpha_size) return false;

  // the multiplication by 1.1 is sometimes necessary, presumably to
  // prevent the upper bound from being too tight:
  if (r_max_min > c_max_min)
    *ub = 1.1 * log(1.0 * (alpha_size - l_r))/r_max_min;
  else
    *ub = 1.1 * log(1.0 * (alpha_size - l_c))/c_max_min;

  return true;
}

bool LambdaCalculator::binary_search(double** matrix, int alpha_size, double lb, double ub, vector<double>& letprob1, vector<double>& letprob2, double* lambda, int maxiter)
{
  double l=0;
  double r=0;
  double l_sum=0;
  double r_sum=0;
  int iter=0;

  while (iter < maxiter && (l>=r || (l_sum < 1.0 && r_sum < 1.0) || (l_sum > 1.0 && r_sum > 1.0)))
    {
      l = lb + (ub - lb)*rand()/RAND_MAX;
      r = lb + (ub - lb)*rand()/RAND_MAX;

      if (!calculate_inv_sum(matrix, alpha_size, l, &l_sum) || !calculate_inv_sum(matrix, alpha_size, r, &r_sum))
        {
          l = 0;
          r = 0;
        }
      iter++;
    }

  if (iter >= maxiter)
    return false;

  while (l_sum != 1.0 && r_sum != 1.0 && (l+r)/2.0 != l && (l+r)/2.0 != r)
    {
      double mid = (l + r)/2.0;
      double mid_sum;
      if (!calculate_inv_sum(matrix, alpha_size, mid, &mid_sum))
        return false;

      if (fabs(mid_sum) >= DBL_MAX)
        return false;

      if ((l_sum < 1.0 && mid_sum >= 1.0) || (l_sum > 1.0 && mid_sum <= 1.0))
        {
          r = mid;
          r_sum = mid_sum;
        }

      else if ((r_sum < 1.0 && mid_sum >= 1.0) || (r_sum > 1.0 && mid_sum <= 1.0))
        {
          l = mid;
          l_sum = mid_sum;
        }

      else
        return false;
    }

  if (fabs(l_sum - 1.0) < fabs(r_sum - 1.0))
    {
      if (check_lambda(matrix, l, alpha_size, letprob1, letprob2))
        {
          *lambda = l;
          return true;
        }
      return false;
    }

  if (check_lambda(matrix, r, alpha_size, letprob1, letprob2))
    {
      *lambda = r;
      return true;
    }
  return false;
}

double LambdaCalculator::calculate_lambda(double** matrix, int alpha_size, vector<double>& letprob1, vector<double>& letprob2, int maxiter, int max_boundary_search_iter, double lb_ratio)
{
  double ub;

  if (!find_ub(matrix, alpha_size, &ub))
    return -1;

  double lb = ub*lb_ratio;
  double lambda = -1;
  int iter = 0;
  bool flag = false;

  while (!flag && iter < maxiter)
    {
      flag = binary_search(matrix, alpha_size, lb, ub, letprob1, letprob2, &lambda, max_boundary_search_iter);
      iter++;
    }

  return lambda;
}

bool LambdaCalculator::check_lambda(double** matrix, double lambda, int alpha_size, vector<double>& letprob1, vector<double>& letprob2)
{
  double **m = makematrix(alpha_size, alpha_size, 0);
  double **y = makematrix(alpha_size, alpha_size, 0);

  for (int i=0; i<alpha_size; i++)
    for (int j=0; j<alpha_size; j++)
      m[i][j] = exp(lambda * matrix[i][j]);

  invert(m, y, alpha_size);

  for (int i=0; i<alpha_size; i++)
    {
      double p = 0;
      for (int j=0;j<alpha_size; j++)
        p += y[i][j];
      if (p < 0 || p > 1)
        {
          letprob2.clear();
          return false;
        }
      letprob2.push_back(roundToFewDigits(p));
    }

  for (int j=0; j<alpha_size; j++)
    {
      double q = 0;
      for (int i=0; i<alpha_size; i++)
        q += y[i][j];
      if (q < 0 || q > 1)
        {
          letprob2.clear();
          letprob1.clear();
          return false;
        }
      letprob1.push_back(roundToFewDigits(q));
    }

  deletematrix(m, alpha_size);
  deletematrix(y, alpha_size);

  return true;
}

void LambdaCalculator::calculate(const const_int_ptr *matrix, int alphSize) {
  assert(alphSize >= 0);
  setBad();

  int maxiter = 1000;
  int max_boundary_search_iter = 100;
  double lb_ratio = 1e-6;

  double** mat = makematrix(alphSize, alphSize, 0);
  for (int i=0; i<alphSize; i++)
    for (int j=0; j<alphSize; j++)
      mat[i][j] = matrix[i][j];
  lambda_ = calculate_lambda(mat, alphSize, letterProbs1_, letterProbs2_, maxiter, max_boundary_search_iter, lb_ratio);
  deletematrix(mat, alphSize);
}
}
