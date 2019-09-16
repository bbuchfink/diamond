// Copyright 2008 Michiaki Hamada
// Modified 2015 Yutaro Konta

// This class calculates the scale factor (lambda), and the letter
// probabilities, that are implicit in a scoring matrix.  The
// calculation might fail for various reasons, putting it into a
// "bad/undefined" state.

// If the score matrix is symmetric, then the two sets of letter
// probabilities should be identical.  With this code, they might
// differ minutely from exact identity.

#ifndef LAMBDA_CALCULATOR_HH
#define LAMBDA_CALCULATOR_HH

#include <vector>

namespace cbrc{

typedef const int *const_int_ptr;

class LambdaCalculator{
 public:
  LambdaCalculator() { setBad(); }

  void calculate(const const_int_ptr *matrix, int alphSize);

  // Put us in the bad/undefined state.
  void setBad();

  // Are we in the bad/undefined state?
  bool isBad() const { return (lambda_ < 0); }

  // The scale factor.  In the bad/undefined state, it is negative.
  double lambda() const { return lambda_; }

  // The probabilities of letters corresponding to matrix rows (1st index).
  // In the bad/undefined state, it is a null pointer.
  const double *letterProbs1() const {return isBad() ? 0 : &letterProbs1_[0];}

  // The probabilities of letters corresponding to matrix columns (2nd index).
  // In the bad/undefined state, it is a null pointer.
  const double *letterProbs2() const {return isBad() ? 0 : &letterProbs2_[0];}

 private:
  double lambda_;
  std::vector<double> letterProbs1_;
  std::vector<double> letterProbs2_;

  bool find_ub(double **matrix, int alpha_size, double *ub);
  bool binary_search(double** matrix, int alpha_size, double lb, double ub, std::vector<double>& letprob1, std::vector<double>& letprob2, double* lambda, int maxiter);
  double calculate_lambda(double** matrix, int alpha_size, std::vector<double>& letprob1, std::vector<double>& letprob2, int maxiter, int max_boundary_search_iter, double lb_ratio);
  bool check_lambda(double** matrix, double lambda, int alpha_size, std::vector<double>& letprob1, std::vector<double>& letprob2);
};

}  // end namespace

#endif
