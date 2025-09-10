#pragma once
#include <immintrin.h>
#include <stdint.h>

#if defined(__FMA__)
#define FMA(a,b,c) _mm256_fmadd_ps((a),(b),(c))
#define FMADD256(a,b,c) _mm256_fmadd_ps((a),(b),(c))
#define FMADD128(a,b,c) _mm_fmadd_ps((a),(b),(c))
#else
#define FMA(a,b,c) _mm256_add_ps(_mm256_mul_ps((a),(b)),(c))
#define FMADD256(a,b,c) _mm256_add_ps(_mm256_mul_ps((a),(b)), (c))
#define FMADD128(a,b,c) _mm_add_ps(_mm_mul_ps((a),(b)), (c))
#endif

static inline __m256 log256_approx_pos(__m256 x) {
    const __m256 min_norm_pos = _mm256_set1_ps(1.17549435e-38f);
    x = _mm256_max_ps(x, min_norm_pos);

    __m256i xi = _mm256_castps_si256(x);
    __m256i ei = _mm256_srli_epi32(xi, 23);
    ei = _mm256_sub_epi32(ei, _mm256_set1_epi32(127));
    __m256  e = _mm256_cvtepi32_ps(ei);

    __m256i mant_i = _mm256_and_si256(xi, _mm256_set1_epi32(0x007FFFFF));
    mant_i = _mm256_or_si256(mant_i, _mm256_set1_epi32(0x3F800000));
    __m256  m = _mm256_castsi256_ps(mant_i);

    __m256 r = _mm256_sub_ps(m, _mm256_set1_ps(1.0f));

    const __m256 c0 = _mm256_set1_ps(-0.006074878f);
    const __m256 c1 = _mm256_set1_ps(0.034418594f);
    const __m256 c2 = _mm256_set1_ps(-0.092313768f);
    const __m256 c3 = _mm256_set1_ps(0.164783493f);
    const __m256 c4 = _mm256_set1_ps(-0.239190713f);
    const __m256 c5 = _mm256_set1_ps(0.331334025f);
    const __m256 c6 = _mm256_set1_ps(-0.499801159f);
    const __m256 c7 = _mm256_set1_ps(0.999991477f);
    const __m256 c8 = _mm256_set1_ps(0.000000091f);

    __m256 y = c0;
    y = FMA(y, r, c1);
    y = FMA(y, r, c2);
    y = FMA(y, r, c3);
    y = FMA(y, r, c4);
    y = FMA(y, r, c5);
    y = FMA(y, r, c6);
    y = FMA(y, r, c7);
    y = FMA(y, r, c8);

    const __m256 ln2 = _mm256_set1_ps(0.69314718056f);
    return FMA(e, ln2, y);
}

static inline __m256 log256_ps_approx(__m256 x) {
    const __m256 min_pos = _mm256_set1_ps(1.0e-38f);
    x = _mm256_max_ps(x, min_pos);

    __m256i xi = _mm256_castps_si256(x);
    __m256i expi = _mm256_srli_epi32(xi, 23);
    expi = _mm256_sub_epi32(expi, _mm256_set1_epi32(127));
    __m256  e = _mm256_cvtepi32_ps(expi);

    const __m256i mant_mask = _mm256_set1_epi32(0x007fffff);
    const __m256i one_bits = _mm256_set1_epi32(0x3f800000);
    __m256i mi = _mm256_or_si256(_mm256_and_si256(xi, mant_mask), one_bits);
    __m256  m = _mm256_castsi256_ps(mi);

    const __m256 sqrt2 = _mm256_set1_ps(1.41421356237f);
    const __m256 half = _mm256_set1_ps(0.5f);
    __m256     m_half = _mm256_mul_ps(m, half);
    __m256     gt_mask = _mm256_cmp_ps(m, sqrt2, _CMP_GT_OQ);
    m = _mm256_blendv_ps(m, m_half, gt_mask);
    e = _mm256_add_ps(e, _mm256_and_ps(gt_mask, _mm256_set1_ps(1.0f)));

    const __m256 one = _mm256_set1_ps(1.0f);
    __m256 f = _mm256_sub_ps(m, one);
    __m256 f2 = _mm256_mul_ps(f, f);

    __m256 sum = f;
    sum = _mm256_sub_ps(sum, _mm256_mul_ps(_mm256_set1_ps(0.5f), f2));
    __m256 t = _mm256_mul_ps(f2, f);
    sum = _mm256_add_ps(sum, _mm256_mul_ps(_mm256_set1_ps(1.0f / 3.0f), t));
    t = _mm256_mul_ps(t, f);
    sum = _mm256_sub_ps(sum, _mm256_mul_ps(_mm256_set1_ps(1.0f / 4.0f), t));
    t = _mm256_mul_ps(t, f);
    sum = _mm256_add_ps(sum, _mm256_mul_ps(_mm256_set1_ps(1.0f / 5.0f), t));
    t = _mm256_mul_ps(t, f);
    sum = _mm256_sub_ps(sum, _mm256_mul_ps(_mm256_set1_ps(1.0f / 6.0f), t));
    t = _mm256_mul_ps(t, f);
    sum = _mm256_add_ps(sum, _mm256_mul_ps(_mm256_set1_ps(1.0f / 7.0f), t));

    const __m256 ln2 = _mm256_set1_ps(0.6931471805599453f);
    return _mm256_add_ps(_mm256_mul_ps(e, ln2), sum);
}

static inline float hsum256_ps(__m256 v) {
    __m128 vlow = _mm256_castps256_ps128(v);
    __m128 vhigh = _mm256_extractf128_ps(v, 1);
    vlow = _mm_add_ps(vlow, vhigh);
    vlow = _mm_hadd_ps(vlow, vlow);
    vlow = _mm_hadd_ps(vlow, vlow);
    return _mm_cvtss_f32(vlow);
}

static inline float hsum128_ps(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);
    __m128 sums = _mm_add_ps(v, shuf);
    shuf = _mm_movehl_ps(shuf, sums);
    sums = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
}

static inline float hsum256_ps_v2(__m256 v) {
    __m128 low = _mm256_castps256_ps128(v);
    __m128 high = _mm256_extractf128_ps(v, 1);
    return hsum128_ps(_mm_add_ps(low, high));
}


static inline float hsum128_ps_v2(__m128 v) {
    __m128 t = _mm_hadd_ps(v, v);
    t = _mm_hadd_ps(t, t);
    return _mm_cvtss_f32(t);
}

static inline float hsum256_ps_v3(__m256 v) {
    __m128 lo = _mm256_castps256_ps128(v);
    __m128 hi = _mm256_extractf128_ps(v, 1);
    return hsum128_ps_v2(_mm_add_ps(lo, hi));
}

static inline float fast_sqrtf_avx(float x) {
    __m256 vx = _mm256_set1_ps(x);
    __m256 r = _mm256_rsqrt_ps(vx);
    __m256 r2 = _mm256_mul_ps(r, r);
    __m256 xr2 = _mm256_mul_ps(vx, r2);
    __m256 three = _mm256_set1_ps(3.0f);
    __m256 half = _mm256_set1_ps(0.5f);
    r = _mm256_mul_ps(r, _mm256_mul_ps(_mm256_sub_ps(three, xr2), half));
    __m256 s = _mm256_mul_ps(vx, r);
    return _mm_cvtss_f32(_mm256_castps256_ps128(s));
}