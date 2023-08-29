// Prefix scan implemenation by Daniel Liu, see: https://github.com/Daniel-Liu-c0deb0t/block-aligner

/*
MIT License

Copyright (c) 2021 Daniel Liu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once
#include "../score_vector_int16.h"

namespace DP { namespace PrefixScan { namespace DISPATCH_ARCH {

#if ARCH_ID == 2 || ARCH_ID == 3

template<int DELTA>
static inline std::pair<::DISPATCH_ARCH::ScoreVector<int16_t, DELTA>, ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA >> prefix_scan_consts(const ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA> gap) {
	__m256i shift1 = _mm256_slli_si256(gap.data_, 2);
	shift1 = _mm256_adds_epi16(shift1, gap.data_);
	__m256i shift2 = _mm256_slli_si256(shift1, 4);
	shift2 = _mm256_adds_epi16(shift2, shift1);
	__m256i shift4 = _mm256_slli_si256(shift2, 8);
	shift4 = _mm256_adds_epi16(shift4, shift2);

	__m256i correct1 = _mm256_srli_si256(_mm256_shufflehi_epi16(shift4, 0xff), 8);
	correct1 = _mm256_permute4x64_epi64(correct1, 5);
	correct1 = _mm256_adds_epi16(correct1, shift4);

	return { ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA>(correct1), ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA>(shift4) };
}

template<int DELTA>
static inline ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA> prefix_scan(
	::DISPATCH_ARCH::ScoreVector<int16_t, DELTA> input,
	const ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA> gap_extend,
	const ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA> gap_cost_lane)
{
	//const __m256i shrt_min = _mm256_set1_epi16(SHRT_MIN);
	__m256i shift1 = _mm256_slli_si256(input.data_, 2);
	//shift1 = _mm256_blend_epi16(shift1, shrt_min, 1);
	shift1 = _mm256_adds_epi16(shift1, gap_extend.data_);
	shift1 = _mm256_max_epi16(input.data_, shift1);
	__m256i shift2 = _mm256_slli_si256(shift1, 4);
	//shift2 = _mm256_blend_epi16(shift2, shrt_min, 3);
	shift2 = _mm256_adds_epi16(shift2, _mm256_slli_epi16(gap_extend.data_, 1));
	shift2 = _mm256_max_epi16(shift1, shift2);
	__m256i shift4 = _mm256_slli_si256(shift2, 8);
	//shift4 = _mm256_blend_epi16(shift4, shrt_min, 15);
	shift4 = _mm256_adds_epi16(shift4, _mm256_slli_epi16(gap_extend.data_, 2));
	shift4 = _mm256_max_epi16(shift2, shift4);
	__m256i correct1 = _mm256_shufflehi_epi16(shift4, 0xff);
	correct1 = _mm256_permute4x64_epi64(correct1, 0x50);
	correct1 = _mm256_adds_epi16(correct1, gap_cost_lane.data_);
	return ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA>(_mm256_max_epi16(shift4, correct1));
}

template<int DELTA>
static inline std::pair<::DISPATCH_ARCH::ScoreVector<int8_t, DELTA>, ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA >> prefix_scan_consts(const ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA> gap) {
	__m256i shift1 = _mm256_slli_si256(gap.data_, 1);
	shift1 = _mm256_adds_epi8(shift1, gap.data_);
	__m256i shift2 = _mm256_slli_si256(shift1, 2);
	shift2 = _mm256_adds_epi8(shift2, shift1);
	__m256i shift4 = _mm256_slli_si256(shift2, 4);
	shift4 = _mm256_adds_epi8(shift4, shift2);
	__m256i shift8 = _mm256_slli_si256(shift4, 8);
	shift8 = _mm256_adds_epi8(shift8, shift4);

	__m256i correct1 = _mm256_srli_si256(_mm256_set1_epi8(_mm256_extract_epi8(shift8, 15)), 8);
	correct1 = _mm256_permute4x64_epi64(correct1, 5);
	correct1 = _mm256_adds_epi8(correct1, shift8);

	return { ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA>(correct1), ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA>(shift8) };
}

template<int DELTA>
static inline ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA> prefix_scan(
	::DISPATCH_ARCH::ScoreVector<int8_t, DELTA> input,
	const ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA> gap_extend,
	const ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA> gap_cost_lane)
{
	const __m256i schar_min = _mm256_set1_epi8(SCHAR_MIN);
	__m256i shift1 = _mm256_slli_si256(input.data_, 1);
	shift1 = _mm256_or_si256(shift1, _mm256_set_epi8('\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', SCHAR_MIN, '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', SCHAR_MIN));
	shift1 = _mm256_adds_epi8(shift1, gap_extend.data_);
	shift1 = _mm256_max_epi8(input.data_, shift1);
	__m256i shift2 = _mm256_slli_si256(shift1, 2);
	//shift2 = _mm256_blendv_epi8(shift2, schar_min, _mm256_set_epi8('\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\xff', '\xff', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\xff', '\xff'));
	shift2 = _mm256_blend_epi16(shift2, schar_min, 1);
	shift2 = _mm256_adds_epi8(shift2, _mm256_set1_epi8(-2));
	shift2 = _mm256_max_epi8(shift1, shift2);
	__m256i shift4 = _mm256_slli_si256(shift2, 4);
	//shift4 = _mm256_blendv_epi8(shift4, schar_min, _mm256_set_epi8('\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\xff', '\xff', '\xff', '\xff', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\xff', '\xff', '\xff', '\xff'));
	shift4 = _mm256_blend_epi16(shift4, schar_min, 3);
	shift4 = _mm256_adds_epi8(shift4, _mm256_set1_epi8(-4));
	shift4 = _mm256_max_epi8(shift2, shift4);
	__m256i shift8 = _mm256_slli_si256(shift4, 8);
	//shift8 = _mm256_blendv_epi8(shift8, schar_min, _mm256_set_epi8('\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff'));
	shift8 = _mm256_blend_epi16(shift8, schar_min, 15);
	shift8 = _mm256_adds_epi8(shift8, _mm256_set1_epi8(-8));
	shift8 = _mm256_max_epi8(shift4, shift8);
	__m256i correct1 = _mm256_set1_epi8(_mm256_extract_epi8(shift8, 15));
	correct1 = _mm256_blendv_epi8(schar_min, correct1, _mm256_set_epi8('\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'));
	correct1 = _mm256_adds_epi8(correct1, gap_cost_lane.data_);
	return ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA>(_mm256_max_epi8(shift8, correct1));
}

#else

template<int DELTA>
static inline std::pair<::DISPATCH_ARCH::ScoreVector<int16_t, DELTA>, ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA >> prefix_scan_consts(const ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA> gap) {
	return { ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA>(), ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA>() };
}

template<int DELTA>
static inline ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA> prefix_scan(::DISPATCH_ARCH::ScoreVector<int16_t, DELTA> input, const ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA> gap_extend, const ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA> gap_cost_lane) {
	return ::DISPATCH_ARCH::ScoreVector<int16_t, DELTA>();
}

template<int DELTA>
static inline std::pair<::DISPATCH_ARCH::ScoreVector<int8_t, DELTA>, ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA >> prefix_scan_consts(const ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA> gap) {
	return { ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA>(), ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA>() };
}

template<int DELTA>
static inline ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA> prefix_scan(::DISPATCH_ARCH::ScoreVector<int8_t, DELTA> input, const ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA> gap_extend, const ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA> gap_cost_lane) {
	return ::DISPATCH_ARCH::ScoreVector<int8_t, DELTA>();
}

#endif

static inline std::pair<int32_t, int32_t> prefix_scan_consts(int32_t gap) {
	return { -1,-1 };
}

static inline int32_t prefix_scan(int32_t input, int32_t, int32_t) {
	return input;
}

}}}