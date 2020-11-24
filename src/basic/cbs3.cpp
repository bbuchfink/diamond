#include "../lib/Eigen/Dense"

using namespace Eigen;

//static void
//ScaledSymmetricProductA(double** W, const Matrix<float, 20, 20>& diagonal, int alphsize)
//{
//    int rowW, colW;   /* iteration indices over the rows and columns of W */
//    int i, j;         /* iteration indices over characters in the alphabet */
//    int m;            /* The number of rows in A; also the size of W */
//
//    m = 2 * alphsize - 1;
//
//    for (rowW = 0; rowW < m; rowW++) {
//        for (colW = 0; colW <= rowW; colW++) {
//            W[rowW][colW] = 0.0;
//        }
//    }
//    for (i = 0; i < alphsize; i++) {
//        for (j = 0; j < alphsize; j++) {
//            double dd;     /* an individual diagonal element */
//
//            dd = diagonal[i * alphsize + j];
//
//            W[j][j] += dd;
//            if (i > 0) {
//                W[i + alphsize - 1][j] += dd;
//                W[i + alphsize - 1][i + alphsize - 1] += dd;
//            }
//        }
//    }
//}