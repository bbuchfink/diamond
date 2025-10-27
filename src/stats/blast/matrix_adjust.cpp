/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

// #define __AVX2__

#include <cmath>
#include <cstring>
#include <algorithm>
#include "math.h"
#include "matrix_adjust.h"

static constexpr int N   = 20;
static constexpr int N2  = N * N;
static constexpr int MA  = 2 * N - 1;
static constexpr int M   = MA + 1;

static inline MatrixFloat**
Nlm_DenseMatrixNew(size_t nrows,
    size_t ncols)
{
    MatrixFloat** mat;

    mat = (MatrixFloat**)calloc(nrows, sizeof(MatrixFloat*));
    if (mat != NULL) {
        mat[0] = (MatrixFloat*)malloc((size_t)nrows *
            (size_t)ncols * sizeof(MatrixFloat));
        if (mat[0] != NULL) {
            for (size_t i = 1; i < nrows; i++) {
                mat[i] = &mat[0][i * ncols];
            }
        }
        else {
            free(mat);
            mat = NULL;
        }
    }
    return mat;
}


static inline MatrixFloat**
Nlm_LtriangMatrixNew(int n)
{
    int i;
    MatrixFloat** L;
    size_t nelts;

    nelts = ((size_t)n * (n + 1)) / 2;

    L = (MatrixFloat**)calloc(n, sizeof(MatrixFloat*));
    if (L != NULL) {
        L[0] = (MatrixFloat*)calloc(nelts, sizeof(MatrixFloat)); // was malloc
        if (L[0] != NULL) {
            for (i = 1; i < n; i++) {
                L[i] = L[i - 1] + i;
            }
        }
        else {
            free(L);
            L = NULL;
        }
    }
    return L;
}

static inline void
Nlm_DenseMatrixFree(MatrixFloat*** mat)
{
    if (*mat != NULL) {
        free((*mat)[0]);
        free(*mat);
    }
    *mat = NULL;
}

#ifdef __AVX2__

static inline void
Nlm_FactorLtriangPosDef(float** A)
{
    for (int i = 0; i < 40; ++i) {
        float* rowi = A[i];
        for (int j = 0; j < i; ++j) {
            const float* rowj = A[j];

            __m256 acc = _mm256_setzero_ps();
            int k = 0;
            int k_end = j & ~7;
            for (; k < k_end; k += 8) {
                __m256 vi = _mm256_loadu_ps(rowi + k);
                __m256 vj = _mm256_loadu_ps(rowj + k);
#if defined(__FMA__) || defined(__FMA4__)
                acc = _mm256_fmadd_ps(vi, vj, acc);
#else
                acc = _mm256_add_ps(acc, _mm256_mul_ps(vi, vj));
#endif
            }
            float dot = hsum256_ps(acc);
            for (; k < j; ++k) dot += rowi[k] * rowj[k];

            float temp = rowi[j] - dot;
            rowi[j] = temp / rowj[j];
        }

        __m256 acc2 = _mm256_setzero_ps();
        int k2 = 0;
        int k2_end = i & ~7;
        for (; k2 < k2_end; k2 += 8) {
            __m256 vi = _mm256_loadu_ps(rowi + k2);
#if defined(__FMA__) || defined(__FMA4__)
            acc2 = _mm256_fmadd_ps(vi, vi, acc2);
#else
            acc2 = _mm256_add_ps(acc2, _mm256_mul_ps(vi, vi));
#endif
        }
        float sumsq = hsum256_ps(acc2);
        for (; k2 < i; ++k2) sumsq += rowi[k2] * rowi[k2];

        float diag = rowi[i] - sumsq;
        if (diag < 0.0f) diag = 0.0f;

        rowi[i] = fast_sqrtf_avx(diag);
    }
}

#else

static inline void
Nlm_FactorLtriangPosDef(MatrixFloat** A)
{
    int i, j, k;
    MatrixFloat temp;
    for (i = 0; i < MA + 1; i++) {
        for (j = 0; j < i; j++) {
            temp = A[i][j];
            for (k = 0; k < j; k++) {
                temp -= A[i][k] * A[j][k];
            }
            A[i][j] = temp / A[j][j];
        }
        temp = A[i][i];
        for (k = 0; k < i; k++) {
            temp -= A[i][k] * A[i][k];
        }
        A[i][i] = sqrt(temp);
    }
}

#endif

#ifdef __AVX2__

static inline void Nlm_SolveLtriangPosDef(float* __restrict x,
    float* const* __restrict L)
{
    for (int i = 0; i < 40; ++i) {
        const float* row = L[i];

        __m256 acc = _mm256_setzero_ps();
        int j = 0;
        for (; j + 7 < i; j += 8) {
            __m256 lv = _mm256_loadu_ps(row + j);
            __m256 xv = _mm256_loadu_ps(x + j);
#if defined(__FMA__)
            acc = _mm256_fmadd_ps(lv, xv, acc);
#else
            acc = _mm256_add_ps(acc, _mm256_mul_ps(lv, xv));
#endif
        }

        float temp = x[i] - hsum256_ps(acc);
        for (; j < i; ++j) temp -= row[j] * x[j];

        x[i] = temp / row[i];
    }

    for (int j = 39; j >= 0; --j) {
        float* row = L[j];

        x[j] /= row[j];

        __m256 xjb = _mm256_set1_ps(x[j]);
        int i = 0;
        for (; i + 7 < j; i += 8) {
            __m256 xi = _mm256_loadu_ps(x + i);
            __m256 lj = _mm256_loadu_ps(row + i);
#if defined(__FMA__)
            xi = _mm256_fnmadd_ps(lj, xjb, xi);
#else
            xi = _mm256_sub_ps(xi, _mm256_mul_ps(lj, xjb));
#endif
            _mm256_storeu_ps(x + i, xi);
        }
        for (; i < j; ++i) x[i] -= row[i] * x[j];
    }
}

#else

static inline void Nlm_SolveLtriangPosDef(MatrixFloat x[],
    MatrixFloat** L)
{
    int i, j;
    MatrixFloat temp;
    for (i = 0; i < MA + 1; i++) {
        temp = x[i];
        for (j = 0; j < i; j++) {
            temp -= L[i][j] * x[j];
        }
        x[i] = temp / L[i][i];
    }
    for (j = MA; j >= 0; j--) {
        x[j] /= L[j][j];
        for (i = 0; i < j; i++) {
            x[i] -= L[j][i] * x[j];
        }
    }
}

#endif

static inline MatrixFloat
Nlm_EuclideanNorm(const MatrixFloat v[], int n)
{
    MatrixFloat sum = 1.0f;
    MatrixFloat scale = 0.0f;
    int i;

    for (i = 0; i < n; i++) {
        if (v[i] != 0.0f) {
            MatrixFloat absvi = fabs(v[i]);
            if (scale < absvi) {
                sum = 1.0f + sum * (scale / absvi) * (scale / absvi);
                scale = absvi;
            }
            else {
                sum += (absvi / scale) * (absvi / scale);
            }
        }
    }
    return scale * sqrt(sum);
}

static inline void Nlm_AddVectors(MatrixFloat y[], int n, MatrixFloat alpha, const MatrixFloat x[])
{
    int i;

    for (i = 0; i < n; i++) {
        y[i] += alpha * x[i];
    }
}

static inline MatrixFloat
Nlm_StepBound(const MatrixFloat x[], int n, const MatrixFloat step_x[], MatrixFloat max)
{
    int i;
    MatrixFloat alpha = max;

    for (i = 0; i < n; i++) {
        MatrixFloat alpha_i;

        alpha_i = -x[i] / step_x[i];
        if (alpha_i >= 0 && alpha_i < alpha) {
            alpha = alpha_i;
        }
    }
    return alpha;
}

#ifdef __AVX2__

static inline void ScaledSymmetricProductA20(float* Wrows[40], const float Dinv[400]) {
    for (int r = 0; r < 39; ++r) {
        int n = r + 1;
        int k = 0;
        const __m256 z = _mm256_setzero_ps();
        for (; k + 8 <= n; k += 8) _mm256_storeu_ps(Wrows[r] + k, z);
        for (; k < n; ++k) Wrows[r][k] = 0.0f;
    }

    alignas(32) float diag_sum[20];
    {
        const __m256 z = _mm256_setzero_ps();
        _mm256_store_ps(diag_sum + 0, z);
        _mm256_store_ps(diag_sum + 8, z);
        _mm_store_ps(diag_sum + 16, _mm_setzero_ps());
    }

    for (int i = 0; i < 20; ++i) {
        const float* row = Dinv + i * 20;

        const __m256 v0 = _mm256_loadu_ps(row + 0);
        const __m256 v1 = _mm256_loadu_ps(row + 8);
        const __m128 v2 = _mm_loadu_ps(row + 16);

        __m256 ds0 = _mm256_load_ps(diag_sum + 0);
        __m256 ds1 = _mm256_load_ps(diag_sum + 8);
        __m128 ds2 = _mm_load_ps(diag_sum + 16);
        ds0 = _mm256_add_ps(ds0, v0);
        ds1 = _mm256_add_ps(ds1, v1);
        ds2 = _mm_add_ps(ds2, v2);
        _mm256_store_ps(diag_sum + 0, ds0);
        _mm256_store_ps(diag_sum + 8, ds1);
        _mm_store_ps(diag_sum + 16, ds2);

        if (i > 0) {
            float* __restrict dst = Wrows[19 + i];

            __m256 d0 = _mm256_loadu_ps(dst + 0);
            __m256 d1 = _mm256_loadu_ps(dst + 8);
            __m128 d2 = _mm_loadu_ps(dst + 16);
            d0 = _mm256_add_ps(d0, v0);
            d1 = _mm256_add_ps(d1, v1);
            d2 = _mm_add_ps(d2, v2);
            _mm256_storeu_ps(dst + 0, d0);
            _mm256_storeu_ps(dst + 8, d1);
            _mm_storeu_ps(dst + 16, d2);

            const float diag_acc = hsum256_ps_v2(v0) + hsum256_ps_v2(v1) + hsum128_ps(v2);
            dst[19 + i] += diag_acc;
        }
    }

    for (int j = 0; j < 20; ++j) {
        Wrows[j][j] += diag_sum[j];
    }
}

#else

static inline void ScaledSymmetricProductA20(MatrixFloat* Wrows[M], const MatrixFloat Dinv[N2]) {
    for (int r = 0; r < MA; ++r) {
        std::fill(Wrows[r], Wrows[r] + (r + 1), 0.0f);
    }
    for (int i = 0; i < N; ++i) {
        const MatrixFloat* row = &Dinv[i * N];
        MatrixFloat diag_acc = 0.0f;
        for (int j = 0; j < N; ++j) {
            const MatrixFloat dd = row[j];
            Wrows[j][j] += dd;
            if (i > 0) {
                Wrows[19 + i][j] += dd;
                diag_acc += dd;
            }
        }
        if (i > 0) {
            Wrows[19 + i][19 + i] += diag_acc;
        }
    }
}

#endif

#ifdef __AVX2__

static inline void MultiplyByA20(float beta, float y[39], float alpha, const float x[400]) {
    if (beta == 0.0f) {
        __m256 z = _mm256_setzero_ps();
        int k = 0;
        for (; k + 8 <= 39; k += 8) _mm256_storeu_ps(y + k, z);
        for (; k < 39; ++k) y[k] = 0.0f;
    }
    else if (beta != 1.0f) {
        __m256 b = _mm256_set1_ps(beta);
        int k = 0;
        for (; k + 8 <= 39; k += 8) {
            __m256 vy = _mm256_loadu_ps(y + k);
            vy = _mm256_mul_ps(vy, b);
            _mm256_storeu_ps(y + k, vy);
        }
        for (; k < 39; ++k) y[k] *= beta;
    }

    const __m256 a8 = _mm256_set1_ps(alpha);
    const __m128 a4 = _mm_set1_ps(alpha);

    for (int i = 0; i < 20; ++i) {
        const float* xi = x + i * 20;

        __m256 sum256 = _mm256_setzero_ps();
        __m128 sum128 = _mm_setzero_ps();

        __m256 x0 = _mm256_loadu_ps(xi + 0);
        __m256 y0 = _mm256_loadu_ps(y + 0);
        sum256 = _mm256_add_ps(sum256, x0);
#if defined(__FMA__)
        __m256 y0n = _mm256_fmadd_ps(x0, a8, y0);
#else
        __m256 v0 = _mm256_mul_ps(x0, a8);
        __m256 y0n = _mm256_add_ps(y0, v0);
#endif
        _mm256_storeu_ps(y + 0, y0n);

        __m256 x1 = _mm256_loadu_ps(xi + 8);
        __m256 y1 = _mm256_loadu_ps(y + 8);
        sum256 = _mm256_add_ps(sum256, x1);
#if defined(__FMA__)
        __m256 y1n = _mm256_fmadd_ps(x1, a8, y1);
#else
        __m256 v1 = _mm256_mul_ps(x1, a8);
        __m256 y1n = _mm256_add_ps(y1, v1);
#endif
        _mm256_storeu_ps(y + 8, y1n);

        __m128 x2 = _mm_loadu_ps(xi + 16);
        __m128 y2 = _mm_loadu_ps(y + 16);
        sum128 = _mm_add_ps(sum128, x2);
#if defined(__FMA__)
        __m128 y2n = _mm_fmadd_ps(x2, a4, y2);
#else
        __m128 v2 = _mm_mul_ps(x2, a4);
        __m128 y2n = _mm_add_ps(y2, v2);
#endif
        _mm_storeu_ps(y + 16, y2n);

        float row_sum_scaled = alpha * (hsum256_ps_v3(sum256) + hsum128_ps_v2(sum128));
        if (i > 0) y[19 + i] += row_sum_scaled;
    }
}

#else

static inline void MultiplyByA20(MatrixFloat beta, MatrixFloat y[MA], MatrixFloat alpha, const MatrixFloat x[N2]) {
    if (beta == 0.0f) {
        for (int i = 0; i < MA; ++i) y[i] = 0.0f;
    } else if (beta != 1.0f) {
        for (int i = 0; i < MA; ++i) y[i] *= beta;
    }
    for (int i = 0; i < N; ++i) {
        const MatrixFloat* xi = &x[i * N];
        MatrixFloat row_sum = 0.0f;
        for (int j = 0; j < N; ++j) {
            const MatrixFloat v = alpha * xi[j];
            y[j] += v;
            row_sum += v;
        }
        if (i > 0) y[19 + i] += row_sum;
    }
}

#endif

#ifdef __AVX2__

static inline void MultiplyByATranspose20(float beta, float y[400], float alpha, const float x[39]) {
    if (beta == 0.0f) {
        const __m256 z = _mm256_setzero_ps();
        for (int k = 0; k < 400; k += 8) {
            _mm256_storeu_ps(y + k, z);
        }
    }
    else if (beta != 1.0f) {
        const __m256 vbeta = _mm256_set1_ps(beta);
        for (int k = 0; k < 400; k += 8) {
            __m256 vy = _mm256_loadu_ps(y + k);
            vy = _mm256_mul_ps(vy, vbeta);
            _mm256_storeu_ps(y + k, vy);
        }
    }

    const __m256 x0 = _mm256_loadu_ps(x + 0);
    const __m256 x1 = _mm256_loadu_ps(x + 8);
    const __m128 x2 = _mm_loadu_ps(x + 16);

    const __m256 valpha256 = _mm256_set1_ps(alpha);
    const __m128 valpha128 = _mm_set1_ps(alpha);

    for (int i = 0; i < 20; ++i) {
        const float add_row_scalar = (i > 0) ? x[19 + i] : 0.0f;
        const __m256 vadd256 = _mm256_set1_ps(add_row_scalar);
        const __m128 vadd128 = _mm_set1_ps(add_row_scalar);

        float* yi = y + i * 20;

        __m256 y0 = _mm256_loadu_ps(yi + 0);
        __m256 y1 = _mm256_loadu_ps(yi + 8);
        __m128 y2 = _mm_loadu_ps(yi + 16);

        const __m256 t0 = _mm256_add_ps(x0, vadd256);
        const __m256 t1 = _mm256_add_ps(x1, vadd256);
        const __m128 t2 = _mm_add_ps(x2, vadd128);

        y0 = FMADD256(valpha256, t0, y0);
        y1 = FMADD256(valpha256, t1, y1);
        y2 = FMADD128(valpha128, t2, y2);

        _mm256_storeu_ps(yi + 0, y0);
        _mm256_storeu_ps(yi + 8, y1);
        _mm_storeu_ps(yi + 16, y2);
    }
}

#else

static inline void MultiplyByATranspose20(MatrixFloat beta, MatrixFloat y[N2], MatrixFloat alpha, const MatrixFloat x[MA]) {
    if (beta == 0.0f) {
        for (int k = 0; k < N2; ++k) y[k] = 0.0f;
    } else if (beta != 1.0f) {
        for (int k = 0; k < N2; ++k) y[k] *= beta;
    }
    for (int i = 0; i < N; ++i) {
        const MatrixFloat add_row = (i > 0) ? x[19 + i] : 0.0f;
        const MatrixFloat* xcol = x;
        MatrixFloat* yi = &y[i * N];
        for (int j = 0; j < N; ++j) {
            yi[j] += alpha * (xcol[j] + add_row);
        }
    }
}

#endif

static inline void ResidualsLinearConstraints20(MatrixFloat rA[MA], const MatrixFloat x[N2],
                                                const MatrixFloat row_sums[N], const MatrixFloat col_sums[N]) {
    for (int j = 0; j < N; ++j) rA[j] = col_sums[j];
    for (int i = 1; i < N; ++i) rA[19 + i] = row_sums[i];

    MultiplyByA20(1.0f, rA, -1.0f, x);
}

#ifdef __AVX2__

static inline void DualResiduals20(float resids_x[400],
    float* const grads[2],
    const float z[40]) {
    const float eta = z[39];
    const __m256 v_eta = _mm256_set1_ps(eta);
    const float* __restrict g0 = grads[0];
    const float* __restrict g1 = grads[1];
    float* __restrict r = resids_x;
    for (int k = 0; k < 400; k += 8) {
        __m256 v_g0 = _mm256_loadu_ps(g0 + k);
        __m256 v_g1 = _mm256_loadu_ps(g1 + k);
        __m256 v_prod = _mm256_mul_ps(v_eta, v_g1);
        __m256 v_res = _mm256_sub_ps(v_prod, v_g0);
        _mm256_storeu_ps(r + k, v_res);
    }
    MultiplyByATranspose20(1.0f, resids_x, 1.0f, z);
}

#else

static inline void DualResiduals20(MatrixFloat resids_x[N2], MatrixFloat* const grads[2],
                                   const MatrixFloat z[M]) {    
    const MatrixFloat eta = z[MA];
    for (int k = 0; k < N2; ++k) resids_x[k] = -grads[0][k] + eta * grads[1][k];    
    MultiplyByATranspose20(1.0f, resids_x, 1.0f, z);
}

#endif

static inline void CalculateResiduals20(MatrixFloat* rnorm,
                                        MatrixFloat resids_x[N2],
                                        MatrixFloat resids_z[M],
                                        const MatrixFloat values[2],
                                        MatrixFloat* const grads[2],
                                        const MatrixFloat row_sums[N],
                                        const MatrixFloat col_sums[N],
                                        const MatrixFloat x[N2],
                                        const MatrixFloat z[M],
                                        MatrixFloat target_re) {
    DualResiduals20(resids_x, grads, z);
    const MatrixFloat norm_resids_x = Nlm_EuclideanNorm(resids_x, N2);

    ResidualsLinearConstraints20(resids_z, x, row_sums, col_sums);
    MatrixFloat norm_resids_z;
    resids_z[MA] = target_re - values[1];
    norm_resids_z = Nlm_EuclideanNorm(resids_z, MA + 1);

    *rnorm = std::sqrt(norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z);
}

#ifdef __AVX2__

static inline void EvaluateReFunctions20(float values[2], float* const grads[2],
    const float x[N2], const float q[N2], const float scores[N2])
{
    float* __restrict g0 = grads[0];
    float* __restrict g1 = grads[1];

    const int n = N2;
    const __m256 one = _mm256_set1_ps(1.0f);
    const __m256 eps = _mm256_set1_ps(1e-30f);

    __m256 acc_f = _mm256_setzero_ps();
    __m256 acc_re = _mm256_setzero_ps();

    int k = 0;
    for (; k + 8 <= n; k += 8) {
        __m256 vx = _mm256_loadu_ps(&x[k]);
        __m256 vq = _mm256_loadu_ps(&q[k]);
        __m256 vs = _mm256_loadu_ps(&scores[k]);

        vq = _mm256_max_ps(vq, eps);
        __m256 ratio = _mm256_div_ps(vx, vq);

        __m256 t = log256_approx_pos(ratio);
        __m256 u = _mm256_add_ps(t, vs);

        acc_f = FMA(vx, t, acc_f);

        _mm256_storeu_ps(&g0[k], _mm256_add_ps(t, one));

        acc_re = FMA(vx, u, acc_re);
        _mm256_storeu_ps(&g1[k], _mm256_add_ps(u, one));
    }

    auto hsum256 = [](__m256 v) -> float {
        __m128 vlow = _mm256_castps256_ps128(v);
        __m128 vhigh = _mm256_extractf128_ps(v, 1);
        __m128 vsum = _mm_add_ps(vlow, vhigh);
        vsum = _mm_hadd_ps(vsum, vsum);
        vsum = _mm_hadd_ps(vsum, vsum);
        float out;
        _mm_store_ss(&out, vsum);
        return out;
        };

    float f_sum = hsum256(acc_f);
    float re_sum = hsum256(acc_re);

    values[0] = f_sum;
    values[1] = re_sum;
}

#else

static inline void EvaluateReFunctions20(MatrixFloat values[2], MatrixFloat* const grads[2],
    const MatrixFloat x[N2], const MatrixFloat q[N2],
    const MatrixFloat scores[N2]) {
    MatrixFloat f = 0.0f, re = 0.0f;

    for (int k = 0; k < N2; ++k) {
        const MatrixFloat t = std::log(x[k] / q[k]);
        f += x[k] * t;
        grads[0][k] = t + 1.0f;

        const MatrixFloat u = t + scores[k];
        re += x[k] * u;
        grads[1][k] = u + 1.0f;
    }
    values[0] = f; values[1] = re;
}

#endif

#ifndef __AVX2__

static inline void ComputeScoresFromProbs20(MatrixFloat scores[N2],
                                            const MatrixFloat target_freqs[N2],
                                            const MatrixFloat row_freqs[N],
                                            const MatrixFloat col_freqs[N]) {
    for (int i = 0; i < N; ++i) {
        const MatrixFloat ri = row_freqs[i];
        const int base = i * N;
        for (int j = 0; j < N; ++j) {
            scores[base + j] = std::log(target_freqs[base + j] / (ri * col_freqs[j]));
        }
    }
}

#else

static inline void ComputeScoresFromProbs20(float scores[N2],
    const float target_freqs[N2],
    const float row_freqs[N],
    const float col_freqs[N]) {
    const __m256 col0 = _mm256_loadu_ps(col_freqs + 0);
    const __m256 col1 = _mm256_loadu_ps(col_freqs + 8);

    const __m256i tail_mask_i = _mm256_set_epi32(0, 0, 0, 0, -1, -1, -1, -1);
    const __m256  tail_mask = _mm256_castsi256_ps(tail_mask_i);

    __m256 col2_raw = _mm256_maskload_ps(col_freqs + 16, tail_mask_i);
    __m256 col2 = _mm256_blendv_ps(_mm256_set1_ps(1.0f), col2_raw, tail_mask);

    for (int i = 0; i < N; ++i) {
        const float ri = row_freqs[i];
        const int base = i * N;

        __m256 r = _mm256_set1_ps(ri);

        __m256 denom0 = _mm256_mul_ps(r, col0);
        __m256 denom1 = _mm256_mul_ps(r, col1);
        __m256 denom2 = _mm256_mul_ps(r, col2);

        __m256 num0 = _mm256_loadu_ps(target_freqs + base + 0);
        __m256 num1 = _mm256_loadu_ps(target_freqs + base + 8);

        __m256 num2_raw = _mm256_maskload_ps(target_freqs + base + 16, tail_mask_i);
        __m256 num2 = _mm256_blendv_ps(_mm256_set1_ps(1.0f), num2_raw, tail_mask);

        __m256 ratio0 = _mm256_div_ps(num0, denom0);
        __m256 ratio1 = _mm256_div_ps(num1, denom1);
        __m256 ratio2 = _mm256_div_ps(num2, denom2);

        __m256 y0 = log256_approx_pos(ratio0);
        __m256 y1 = log256_approx_pos(ratio1);
        __m256 y2 = log256_approx_pos(ratio2);

        _mm256_storeu_ps(scores + base + 0, y0);
        _mm256_storeu_ps(scores + base + 8, y1);
        _mm256_maskstore_ps(scores + base + 16, tail_mask_i, y2);
    }
}

#endif

struct NewtonSys20 {
    MatrixFloat Wbuf[(M * (M + 1)) / 2];
    MatrixFloat* Wrows[M];
    MatrixFloat Dinv[N2];
    MatrixFloat grad_re[N2];

    NewtonSys20() {
        size_t off = 0;
        for (int i = 0; i < M; ++i) {
            Wrows[i] = &Wbuf[off];
            off += (i + 1);
        }
    }
};

static inline void FactorNewton20(NewtonSys20& sys,
    const MatrixFloat x[N2],
    const MatrixFloat z[M],
    MatrixFloat* const grads[2],
    MatrixFloat workspace[N2]) {
    const MatrixFloat eta = z[MA];
    const MatrixFloat s = 1.0f / (1.0f - eta);
    for (int k = 0; k < N2; ++k) sys.Dinv[k] = x[k] * s;

    ScaledSymmetricProductA20(sys.Wrows, sys.Dinv);

    std::memcpy(sys.grad_re, grads[1], N2 * sizeof(MatrixFloat));

    MatrixFloat* lastRow = sys.Wrows[MA];
    lastRow[MA] = 0.0f;

    for (int k = 0; k < N2; ++k) workspace[k] = sys.Dinv[k] * sys.grad_re[k];
    for (int k = 0; k < N2; ++k) lastRow[MA] += sys.grad_re[k] * workspace[k];

    MultiplyByA20(0.0f, lastRow, 1.0f, workspace);
    Nlm_FactorLtriangPosDef(sys.Wrows);
}

static inline void SolveNewton20(MatrixFloat step_x[N2], MatrixFloat step_z[M],
                                 const NewtonSys20& sys,
                                 MatrixFloat workspace[N2]) {
    for (int k = 0; k < N2; ++k) workspace[k] = step_x[k] * sys.Dinv[k];
    MultiplyByA20(1.0f, step_z, -1.0f, workspace);
    for (int k = 0; k < N2; ++k) step_z[MA] -= sys.grad_re[k] * workspace[k];
    
    Nlm_SolveLtriangPosDef(step_z, const_cast<MatrixFloat**>(sys.Wrows));

    for (int k = 0; k < N2; ++k) step_x[k] += sys.grad_re[k] * step_z[MA];
    MultiplyByATranspose20(1.0f, step_x, 1.0f, step_z);
    for (int k = 0; k < N2; ++k) step_x[k] *= sys.Dinv[k];
}

int New_OptimizeTargetFrequencies(MatrixFloat x[],
                                    int alphsize,
                                    int* iterations,
                                    const MatrixFloat q[],
                                    const MatrixFloat row_sums[],
                                    const MatrixFloat col_sums[],
                                    MatrixFloat relative_entropy,
                                    MatrixFloat tol,
                                    int maxits)
{
    if (alphsize != N || !x || !q || !row_sums || !col_sums || !iterations) {
        if (iterations) *iterations = 0;
        return -1;
    }
#ifdef __AVX2__
	tol = std::max(tol, 0.00001f);
#endif
    MatrixFloat values[2] = {0.0f, 0.0f};
    MatrixFloat grads0[N2], grads1[N2];
    MatrixFloat* grads[2] = {grads0, grads1};

    NewtonSys20 sys;
    MatrixFloat z[M] = {0.0f};
    MatrixFloat resids_x[N2];
    MatrixFloat resids_z[M];
    MatrixFloat old_scores[N2];
    MatrixFloat workspace[N2];

    ComputeScoresFromProbs20(old_scores, q, row_sums, col_sums);

    std::memcpy(x, q, N2 * sizeof(MatrixFloat));

    int its = 0;
    MatrixFloat rnorm = 0.0f;
    
    while (its <= maxits) {
        EvaluateReFunctions20(values, grads, x, q, old_scores);
        CalculateResiduals20(&rnorm, resids_x, resids_z, values, grads,
                             row_sums, col_sums, x, z, relative_entropy);

        if (!(rnorm > tol)) break;

        if (++its <= maxits) {
            FactorNewton20(sys, x, z, grads, workspace);
            SolveNewton20(resids_x, resids_z, sys, workspace);

            MatrixFloat alpha = Nlm_StepBound(x, N2, resids_x, 1.0f / .95f);
            alpha *= 0.95f;

            Nlm_AddVectors(x, N2, alpha, resids_x);
            Nlm_AddVectors(z, MA + 1, alpha, resids_z);
        }
    }

    int converged = 0;
    if (its <= maxits && rnorm <= tol) {
        if (z[MA] < 1.0f) converged = 1;
    }
    *iterations = its;
    return converged ? 0 : 1;
}