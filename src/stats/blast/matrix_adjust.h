using MatrixFloat = float;

int
New_OptimizeTargetFrequencies(MatrixFloat x[],
    int alphsize,
    int* iterations,
    const MatrixFloat q[],
    const MatrixFloat row_sums[],
    const MatrixFloat col_sums[],
    MatrixFloat relative_entropy,
    MatrixFloat tol,
    int maxits);

namespace Stats {

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
    int maxits);

}