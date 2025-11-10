// LambdaCalculator.cpp (single-file, C++17)
// Computes lambda and letter probabilities for a score matrix S,
// solving sum( inv( exp(lambda * S) ) ) = 1.
//
// Public API:
//   LambdaCalculator::Result r = LambdaCalculator().compute(scores);
//
// Where `scores` is an N x N matrix of ints (alphabet size N).
// On success: r.ok == true, r.lambda set, and r.leftProbs / r.rightProbs set.
//
// No third-party dependencies; uses an in-house LU (partial pivoting).

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

class LambdaCalculator {
public:
    struct Result {
        bool ok = false;
        double lambda = -1.0;
        std::vector<double> leftProbs;   // column sums of inv(M)
        std::vector<double> rightProbs;  // row sums of inv(M)
        std::string reason;
    };

    // Main entry: scores as N x N ints.
    Result compute(const std::vector<std::vector<int>>& scores,
        int maxOuterIters = 1000,
        int maxBracketTries = 100,
        double lbRatio = 1e-6) {
        Result out;
        const int n = static_cast<int>(scores.size());
        if (n == 0) { out.reason = "Empty matrix."; return out; }
        for (const auto& row : scores)
            if (static_cast<int>(row.size()) != n) {
                out.reason = "Matrix must be square."; return out;
            }

        // 1) Find an upper bound for lambda.
        double ub = 0.0;
        if (!findUpperBound(scores, ub)) {
            out.reason = "Failed to find a valid upper bound (score matrix violates sign conditions).";
            return out;
        }
        const double lb = lbRatio * ub;

        // 2) Try randomized bracketing + bisection a few times.
        std::mt19937_64 rng(0xC001D00D);
        std::uniform_real_distribution<double> U(lb, ub);

        for (int attempt = 0; attempt < maxOuterIters && !out.ok; ++attempt) {
            // Randomly search for l,r with invSum(l) and invSum(r) on opposite sides of 1.
            double L = 0.0, R = 0.0;
            double fL = 0.0, fR = 0.0; // f(λ) := sum(inv(M(λ))) - 1
            bool bracketed = false;

            for (int t = 0; t < maxBracketTries && !bracketed; ++t) {
                L = U(rng);
                R = U(rng);
                if (L > R) std::swap(L, R);

                if (!invSum(scores, L, fL)) continue;
                if (!invSum(scores, R, fR)) continue;

                if ((fL <= 0.0 && fR >= 0.0) || (fL >= 0.0 && fR <= 0.0)) {
                    bracketed = true;
                }
            }

            if (!bracketed) continue;

            // Bisection on [L, R] until convergence.
            const double relTol = 1e-12;
            const double absTolF = 1e-12;
            int safeGuard = 200; // bisection steps
            while (safeGuard-- > 0) {
                const double mid = 0.5 * (L + R);
                double fM = 0.0;
                if (!invSum(scores, mid, fM)) { out.reason = "Singular/unstable matrix during bisection."; break; }

                if (std::abs(fM) <= absTolF) {
                    // Found a λ that makes sum == 1 (within tolerance)
                    if (finalize(scores, mid, out)) return out;
                    break;
                }

                // Stop if interval is too tight to move further.
                if (!(mid > L && mid < R)) break;
                if ((fM < 0.0 && fL < 0.0) || (fM > 0.0 && fL > 0.0)) {
                    L = mid; fL = fM;
                }
                else {
                    R = mid; fR = fM;
                }

                if (std::abs(R - L) <= std::max(1.0, std::abs(mid)) * relTol) {
                    // Choose the better endpoint.
                    const double lambdaCand = (std::abs(fL) < std::abs(fR)) ? L : R;
                    if (finalize(scores, lambdaCand, out)) return out;
                    break;
                }
            }

            // If bisection loop exits without success, try another random bracket.
        }

        if (!out.ok && out.reason.empty())
            out.reason = "Failed to bracket and solve for lambda.";
        return out;
    }

private:
    // Small row-major matrix helper on a flat vector<double>.
    struct Mat {
        int n = 0;
        std::vector<double> a; // size n*n
        Mat() = default;
        explicit Mat(int n_) : n(n_), a(static_cast<size_t>(n_)* n_, 0.0) {}
        inline double& at(int i, int j) { return a[static_cast<size_t>(i) * n + j]; }
        inline const double& at(int i, int j) const { return a[static_cast<size_t>(i) * n + j]; }
    };

    // Tidy a floating value by formatting to a reasonable precision then parsing.
    static double tidy(double x) {
        std::ostringstream oss;
        oss.setf(std::ios::fmtflags(0), std::ios::floatfield);
        oss << std::setprecision(8) << x; // a few significant digits to collapse tiny noise
        std::istringstream iss(oss.str());
        double y = x;
        iss >> y;
        return y;
    }

    // Build M(lambda) = exp(lambda * S), S integer matrix.
    static Mat buildExpMatrix(const std::vector<std::vector<int>>& S, double lambda) {
        const int n = static_cast<int>(S.size());
        Mat M(n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                M.at(i, j) = std::exp(lambda * static_cast<double>(S[i][j]));
        return M;
    }

    // LU (Doolittle) with partial pivoting, in-place on A.
    // A becomes L+U (unit diag L). Returns false if singular.
    static bool luDecompose(Mat& A, std::vector<int>& piv, double eps = 1e-12) {
        const int n = A.n;
        piv.resize(n);
        for (int i = 0; i < n; ++i) piv[i] = i;

        for (int k = 0; k < n; ++k) {
            // Pivot search on column k, rows k..n-1
            int p = k;
            double maxAbs = std::abs(A.at(k, k));
            for (int i = k + 1; i < n; ++i) {
                double v = std::abs(A.at(i, k));
                if (v > maxAbs) { maxAbs = v; p = i; }
            }
            if (maxAbs < eps) return false; // singular/ill-conditioned for our purposes

            if (p != k) {
                // swap rows
                for (int j = 0; j < n; ++j)
                    std::swap(A.at(k, j), A.at(p, j));
                std::swap(piv[k], piv[p]);
            }

            // Elimination
            for (int i = k + 1; i < n; ++i) {
                A.at(i, k) /= A.at(k, k);
                const double lik = A.at(i, k);
                for (int j = k + 1; j < n; ++j) {
                    A.at(i, j) -= lik * A.at(k, j);
                }
            }
        }
        return true;
    }

    // Solve A x = b using LU (A is LU-packed, piv is permutation).
    static void luSolve(const Mat& LU, const std::vector<int>& piv,
        const std::vector<double>& b,
        std::vector<double>& x) {
        const int n = LU.n;
        x.assign(n, 0.0);

        // Apply permutation to b: Pb
        std::vector<double> y(n);
        for (int i = 0; i < n; ++i) y[i] = b[piv[i]];

        // Forward solve L z = Pb (L unit-diagonal)
        for (int i = 0; i < n; ++i) {
            double sum = y[i];
            for (int j = 0; j < i; ++j) sum -= LU.at(i, j) * y[j];
            y[i] = sum;
        }

        // Backward solve U x = z
        for (int i = n - 1; i >= 0; --i) {
            double sum = y[i];
            for (int j = i + 1; j < n; ++j) sum -= LU.at(i, j) * x[j];
            x[i] = sum / LU.at(i, i);
        }
    }

    // Invert a matrix via LU. Returns false if singular.
    static bool invert(const Mat& A, Mat& Ainv) {
        Mat LU = A;
        std::vector<int> piv;
        if (!luDecompose(LU, piv)) return false;

        const int n = A.n;
        Ainv = Mat(n);
        std::vector<double> e(n), col(n);

        for (int j = 0; j < n; ++j) {
            std::fill(e.begin(), e.end(), 0.0);
            e[j] = 1.0;
            luSolve(LU, piv, e, col);
            for (int i = 0; i < n; ++i) Ainv.at(i, j) = col[i];
        }
        return true;
    }

    // Compute f(λ) = sum(inv(M(λ))) - 1. Returns false if inversion fails.
    static bool invSum(const std::vector<std::vector<int>>& S, double lambda, double& f) {
        Mat M = buildExpMatrix(S, lambda);
        Mat Minv;
        if (!invert(M, Minv)) return false;

        long double acc = 0.0L;
        const int n = M.n;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                acc += static_cast<long double>(Minv.at(i, j));

        const long double val = acc - 1.0L;
        if (!std::isfinite(static_cast<double>(val))) return false;
        f = static_cast<double>(val);
        return true;
    }

    // Once a lambda is found (or a candidate chosen), compute probs and validate.
    bool finalize(const std::vector<std::vector<int>>& S, double lambda, Result& out) const {
        Mat M = buildExpMatrix(S, lambda);
        Mat Minv;
        if (!invert(M, Minv)) { out.reason = "Matrix inversion failed at finalization."; return false; }

        const int n = M.n;
        std::vector<double> rowSums(n, 0.0), colSums(n, 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                rowSums[i] += Minv.at(i, j);
                colSums[j] += Minv.at(i, j);
            }

        // Validate probabilities in [0,1].
        for (double v : rowSums)
            if (!(v >= 0.0 && v <= 1.0)) { out.reason = "Row probability outside [0,1]."; return false; }
        for (double v : colSums)
            if (!(v >= 0.0 && v <= 1.0)) { out.reason = "Column probability outside [0,1]."; return false; }

        // Light rounding to collapse fp jitter (e.g., to make 0.25 exactly).
        for (double& v : rowSums) v = tidy(v);
        for (double& v : colSums) v = tidy(v);

        out.ok = true;
        out.lambda = lambda;
        out.rightProbs = std::move(rowSums);
        out.leftProbs = std::move(colSums);
        out.reason.clear();
        return true;
    }

    // Upper bound heuristic:
    // For each nonzero row i, require min<0<max and track rMaxMin := min over rows of (row max).
    // Same for columns -> cMaxMin. Count all-zero rows/cols (lr, lc).
    // Use 1.1 * log(N - lr) / rMaxMin or with columns, picking the tighter one.
    static bool findUpperBound(const std::vector<std::vector<int>>& S, double& ub) {
        const int n = static_cast<int>(S.size());
        const double INF = std::numeric_limits<double>::infinity();

        double rMaxMin = INF;
        double cMaxMin = INF;
        int lr = 0, lc = 0;

        // Rows
        for (int i = 0; i < n; ++i) {
            int rmax = std::numeric_limits<int>::min();
            int rmin = std::numeric_limits<int>::max();
            for (int j = 0; j < n; ++j) {
                rmax = std::max(rmax, S[i][j]);
                rmin = std::min(rmin, S[i][j]);
            }
            if (rmax == 0 && rmin == 0) { ++lr; continue; }
            if (rmax <= 0 || rmin >= 0) return false; // no sign change
            rMaxMin = std::min(rMaxMin, static_cast<double>(rmax));
        }

        // Cols
        for (int j = 0; j < n; ++j) {
            int cmax = std::numeric_limits<int>::min();
            int cmin = std::numeric_limits<int>::max();
            for (int i = 0; i < n; ++i) {
                cmax = std::max(cmax, S[i][j]);
                cmin = std::min(cmin, S[i][j]);
            }
            if (cmax == 0 && cmin == 0) { ++lc; continue; }
            if (cmax <= 0 || cmin >= 0) return false;
            cMaxMin = std::min(cMaxMin, static_cast<double>(cmax));
        }

        if (lr == n) return false; // all-zero matrix

        // Slight expansion to avoid being too tight, mirroring the intent of 1.1 factor.
        const double N_eff_rows = static_cast<double>(n - lr);
        const double N_eff_cols = static_cast<double>(n - lc);

        if (rMaxMin > cMaxMin) {
            ub = 1.1 * std::log(N_eff_rows) / rMaxMin;
        }
        else {
            ub = 1.1 * std::log(N_eff_cols) / cMaxMin;
        }
        return std::isfinite(ub) && ub > 0.0;
    }
};

// --------------------------
// Example usage / quick test
// --------------------------
#ifdef LAMBDA_CALCULATOR_DEMO
int main() {
    // Toy example: a simple DNA-like 4x4 matrix with positive match and negative mismatch.
    // This is just to sanity-check execution; real inputs come from scoring schemes.
    std::vector<std::vector<int>> S = {
      {  2, -1, -1, -1 },
      { -1,  2, -1, -1 },
      { -1, -1,  2, -1 },
      { -1, -1, -1,  2 }
    };

    LambdaCalculator::Result res = LambdaCalculator().compute(S);
    if (!res.ok) {
        std::cerr << "Failed: " << res.reason << "\n";
        return 1;
    }
    std::cout << "lambda = " << std::setprecision(12) << res.lambda << "\n";
    std::cout << "left probs: ";
    for (double v : res.leftProbs) std::cout << v << " ";
    std::cout << "\nright probs: ";
    for (double v : res.rightProbs) std::cout << v << " ";
    std::cout << "\n";
    return 0;
}
#endif