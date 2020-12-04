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