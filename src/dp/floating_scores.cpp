// SPDX-License-Identifier: GPL-3.0-or-later
#include "floating_scores.h"
#include "util/simd.h"
#include "util/simd/dispatch.h"
#include <array>
#include <algorithm>
#include <limits>
#include <numeric>

namespace DP { namespace DISPATCH_ARCH {

// Double accumulation keeps every bit of the existing float corrections and
// matches the scalar reference, without quantizing to fixed-point scores.
struct DoubleVector {
#if ARCH_AVX2_KERNELS
	static constexpr size_t LANES = 4;
	__m256d value;
	explicit DoubleVector(double x) : value(_mm256_set1_pd(x)) {}
	explicit DoubleVector(const double* x) : value(_mm256_loadu_pd(x)) {}
	explicit DoubleVector(__m256d x) : value(x) {}
	void store(double* x) const { _mm256_storeu_pd(x, value); }
	DoubleVector operator+(DoubleVector x) const { return DoubleVector(_mm256_add_pd(value, x.value)); }
	DoubleVector operator-(DoubleVector x) const { return DoubleVector(_mm256_sub_pd(value, x.value)); }
	DoubleVector max(DoubleVector x) const { return DoubleVector(_mm256_max_pd(value, x.value)); }
#elif defined(__SSE2__)
	static constexpr size_t LANES = 2;
	__m128d value;
	explicit DoubleVector(double x) : value(_mm_set1_pd(x)) {}
	explicit DoubleVector(const double* x) : value(_mm_loadu_pd(x)) {}
	explicit DoubleVector(__m128d x) : value(x) {}
	void store(double* x) const { _mm_storeu_pd(x, value); }
	DoubleVector operator+(DoubleVector x) const { return DoubleVector(_mm_add_pd(value, x.value)); }
	DoubleVector operator-(DoubleVector x) const { return DoubleVector(_mm_sub_pd(value, x.value)); }
	DoubleVector max(DoubleVector x) const { return DoubleVector(_mm_max_pd(value, x.value)); }
#elif defined(__aarch64__)
	static constexpr size_t LANES = 2;
	float64x2_t value;
	explicit DoubleVector(double x) : value(vdupq_n_f64(x)) {}
	explicit DoubleVector(const double* x) : value(vld1q_f64(x)) {}
	explicit DoubleVector(float64x2_t x) : value(x) {}
	void store(double* x) const { vst1q_f64(x, value); }
	DoubleVector operator+(DoubleVector x) const { return DoubleVector(vaddq_f64(value, x.value)); }
	DoubleVector operator-(DoubleVector x) const { return DoubleVector(vsubq_f64(value, x.value)); }
	DoubleVector max(DoubleVector x) const { return DoubleVector(vmaxq_f64(value, x.value)); }
#else
	static constexpr size_t LANES = 1;
	double value;
	explicit DoubleVector(double x) : value(x) {}
	explicit DoubleVector(const double* x) : value(*x) {}
	void store(double* x) const { *x = value; }
	DoubleVector operator+(DoubleVector x) const { return DoubleVector(value + x.value); }
	DoubleVector operator-(DoubleVector x) const { return DoubleVector(value - x.value); }
	DoubleVector max(DoubleVector x) const { return DoubleVector(std::max(value, x.value)); }
#endif
};

std::vector<double> floating_scores(Sequence query, const std::vector<float>& bias,
	const std::vector<Sequence>& targets, FloatingScoring scoring) {
	constexpr size_t lanes = DoubleVector::LANES;
	const double negative = -std::numeric_limits<double>::infinity();
	std::vector<double> result(targets.size(), 0.0);
	std::vector<size_t> order(targets.size());
	std::iota(order.begin(), order.end(), 0);
	// Similar lengths reduce work in lanes whose target has already ended.
	std::stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
		return targets[i].length() < targets[j].length();
		});
	const DoubleVector zero(0.0), gap_open(double(scoring.gap_open)), gap_extend(double(scoring.gap_extend));
	for (size_t begin = 0; begin < order.size(); begin += lanes) {
		const size_t count = std::min(lanes, order.size() - begin);
		const size_t length = targets[order[begin + count - 1]].length();
		const size_t width = (length + 1) * lanes;
		std::vector<double> previous(width, 0.0), current(width, 0.0), vertical(width, negative);
		std::array<std::vector<double>, 32> profiles;
		DoubleVector best(0.0);
		for (size_t i = 0; i < query.length(); ++i) {
			auto& profile = profiles[query[i]];
			if (profile.empty()) {
				profile.assign(width, negative);
				for (size_t j = 1; j <= length; ++j)
					for (size_t lane = 0; lane < count; ++lane) {
						const Sequence target = targets[order[begin + lane]];
						if (j <= target.length())
							profile[j * lanes + lane] = scoring.matrix[size_t(query[i]) * 32 + target[j - 1]];
					}
			}
			const DoubleVector correction{ double(bias[i]) };
			DoubleVector horizontal(negative);
			zero.store(current.data());
			for (size_t j = 1; j <= length; ++j) {
				const size_t p = j * lanes;
				horizontal = (DoubleVector(current.data() + p - lanes) - gap_open - gap_extend).max(horizontal - gap_extend);
				const DoubleVector v = (DoubleVector(previous.data() + p) - gap_open - gap_extend)
					.max(DoubleVector(vertical.data() + p) - gap_extend);
				const DoubleVector diagonal = DoubleVector(previous.data() + p - lanes)
					+ (DoubleVector(profile.data() + p) + correction);
				const DoubleVector score = zero.max(diagonal).max(horizontal).max(v);
				score.store(current.data() + p);
				v.store(vertical.data() + p);
				best = best.max(score);
			}
			previous.swap(current);
		}
		double scores[lanes];
		best.store(scores);
		for (size_t lane = 0; lane < count; ++lane) result[order[begin + lane]] = scores[lane];
	}
	return result;
}

}
DISPATCH_4(std::vector<double>, floating_scores, Sequence, query, const std::vector<float>&, bias, const std::vector<Sequence>&, targets, FloatingScoring, scoring)
}
