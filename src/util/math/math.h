/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <cstdint>
#include <cstring>

namespace Math {

inline double from_bits(std::uint64_t b) {
    double d;
    std::memcpy(&d, &b, sizeof d);
    return d;
}

inline std::uint64_t to_bits(double d) {
    std::uint64_t b;
    std::memcpy(&b, &d, sizeof d);
    return b;
}

inline double quiet_nan() { return from_bits(0x7ff8000000000000ULL); }

inline double exp(double x) {
    constexpr double LN2_HI  = 6.93147180369123816490e-01;
    constexpr double LN2_LO  = 1.90821492927058770002e-10;
    constexpr double INV_LN2 = 1.44269504088896338700e+00;    
    constexpr double RND_SHIFT = 6755399441055744.0;

    const std::uint64_t ux  = to_bits(x);
    const std::uint64_t ax  = ux & 0x7fffffffffffffffULL;

    if (ax > 0x7ff0000000000000ULL) return quiet_nan();
    if (ax == 0x7ff0000000000000ULL)
        return (ux >> 63) ? 0.0 : x;
    if (x >= 709.78271289338409) return from_bits(0x7ff0000000000000ULL);
    if (x <= -745.13321910194122) return 0.0;

    double t  = x * INV_LN2 + RND_SHIFT;
    double kd = t - RND_SHIFT;
    int    k  = static_cast<int>(kd);
    double r  = x - kd * LN2_HI;
    r         = r - kd * LN2_LO;

    const double c2  = 5.00000000000000000000e-01;
    const double c3  = 1.66666666666666666667e-01;
    const double c4  = 4.16666666666666666667e-02;
    const double c5  = 8.33333333333333333333e-03;
    const double c6  = 1.38888888888888888889e-03;
    const double c7  = 1.98412698412698412698e-04;
    const double c8  = 2.48015873015873015873e-05;
    const double c9  = 2.75573192239858906526e-06;
    const double c10 = 2.75573192239858906526e-07;
    const double c11 = 2.50521083854417187751e-08;
    const double c12 = 2.08767569878680989792e-09;
    const double c13 = 1.60590438368216145994e-10;

    double p = c13;
    p = c12 + r * p;
    p = c11 + r * p;
    p = c10 + r * p;
    p = c9  + r * p;
    p = c8  + r * p;
    p = c7  + r * p;
    p = c6  + r * p;
    p = c5  + r * p;
    p = c4  + r * p;
    p = c3  + r * p;
    p = c2  + r * p;
    double er = 1.0 + r + (r * r) * p;

    if (k >= -1021 && k <= 1023) {
        return er * from_bits(static_cast<std::uint64_t>(1023 + k) << 52);
    }
    if (k > 1023) {
        er *= from_bits(static_cast<std::uint64_t>(1023 + 1023) << 52);
        k  -= 1023;
        return er * from_bits(static_cast<std::uint64_t>(1023 + k) << 52);
    }
    er *= from_bits(static_cast<std::uint64_t>(1023 - 1000) << 52);
    k  += 1000;
    return er * from_bits(static_cast<std::uint64_t>(1023 + k) << 52);
}

inline double log(double x) {
    constexpr double LN2_HI = 6.93147180369123816490e-01;
    constexpr double LN2_LO = 1.90821492927058770002e-10;
    constexpr double SQRT2  = 1.41421356237309514547e+00;

    const double Lg1 = 6.666666666666735130e-01;
    const double Lg2 = 3.999999999940941908e-01;
    const double Lg3 = 2.857142874366239149e-01;
    const double Lg4 = 2.222219843214978396e-01;
    const double Lg5 = 1.818357216161805012e-01;
    const double Lg6 = 1.531383769920937332e-01;
    const double Lg7 = 1.479819860511658591e-01;

    std::uint64_t u = to_bits(x);

    if (u >> 63) {
        if ((u << 1) == 0) return from_bits(0xfff0000000000000ULL);
        return quiet_nan();
    }
    if (u == 0) return from_bits(0xfff0000000000000ULL);
    if (u >= 0x7ff0000000000000ULL) {
        if (u == 0x7ff0000000000000ULL) return x;
        return quiet_nan();
    }

    int e = 0;
    if (u < 0x0010000000000000ULL) {
        x *= 18014398509481984.0;
        u  = to_bits(x);
        e -= 54;
    }

    e += static_cast<int>(u >> 52) - 1023;
    u  = (u & 0x000fffffffffffffULL) | 0x3ff0000000000000ULL;
    // Halves a mantissa greater than sqrt(2) by decrementing the exponent, without a hardly predictable branch.
    const std::uint64_t high = u > to_bits(SQRT2) ? 1 : 0;
    u -= high << 52;
    e += static_cast<int>(high);
    double m = from_bits(u);

    double f = m - 1.0;
    double s = f / (2.0 + f);
    double z = s * s;
    double w = z * z;
    double t1 = w * (Lg2 + w * (Lg4 + w * Lg6));
    double t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
    double R  = t2 + t1;
    double hfsq = 0.5 * f * f;

    double dk = static_cast<double>(e);
    return dk * LN2_HI - ((hfsq - (s * (hfsq + R) + dk * LN2_LO)) - f);
}

}