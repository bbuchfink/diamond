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

#include "cbs.h"
#include "basic/config.h"
#include "score_matrix.h"
#include "util/log_stream.h"

using std::array;
using std::runtime_error;
using std::vector;

namespace Stats {

CBS comp_based_stats(0, -1.0, -1.0, -1.0);

CBS::CBS(unsigned code, double query_match_distance_threshold, double length_ratio_threshold, double angle):
    query_match_distance_threshold(-1.0),
    length_ratio_threshold(-1.0),
    angle(50.0)
{
    /*switch (code) {
    case COMP_BASED_STATS_AND_MATRIX_ADJUST:
        this->angle = 70.0;
        this->query_match_distance_threshold = 0.16;
        this->length_ratio_threshold = 3.0;
    default:
        ;
    }*/
    if (angle != -1.0)
        this->angle = angle;
    if (query_match_distance_threshold != 1.0)
        this->query_match_distance_threshold = query_match_distance_threshold;
    if (length_ratio_threshold != -1.0)
        this->length_ratio_threshold = length_ratio_threshold;
}

Composition composition(const Sequence& s) {
    Composition r;
    r.fill(0.0);
    int n = 0;
    for (Loc i = 0; i < s.length(); ++i) {
        int l = s[i];
        if (l < 20) {
            ++r[l];
            ++n;
        }
    }
    if (n == 0)
        return r;
    for (int i = 0; i < 20; ++i)
        r[i] /= n;
    return r;
}

int count_true_aa(const Sequence& s) {
    int n = 0;
    for (Loc i = 0; i < s.length(); ++i)
        if ((size_t)s[i] < TRUE_AA)
            ++n;
    return n;
}

bool use_seg_masking(const Sequence& a, const Sequence& b) {
    //if (config.comp_based_stats != CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST || a.length() != b.length())
    return true;
    Loc n = 0;
    for (Loc i = 0; i < a.length(); ++i)
        if (a[i] == b[i])
            ++n;
    return n != a.length();
}

int TargetMatrix::score_width() const {
    return (score_max > SCHAR_MAX || score_min < SCHAR_MIN) ? 1 : 0;
}

EMatrixAdjustRule adjust_matrix(const Composition& query_comp, int query_len, unsigned cbs, const Sequence& target) {
    if (!CBS::matrix_adjust(cbs) || target.length() == 0 || query_len == 0)
        return eDontAdjustMatrix;

    //Masking::get()(target_seq.data(), target_seq.size(), Masking::Algo::SEG);

    const auto c = composition(target);
    if (CBS::conditioned(cbs)) {
        const EMatrixAdjustRule rule = s_TestToApplyREAdjustmentConditional(query_len, (int)target.length(), query_comp.data(), c.data(), score_matrix.background_freqs());
        if (cbs == CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST)
            return rule;
        else
			return rule == eUserSpecifiedRelEntropy ? eUserSpecifiedRelEntropy : eDontAdjustMatrix;
    }
    else
        return eUserSpecifiedRelEntropy;
}

TargetMatrix::TargetMatrix(const Composition& query_comp, int query_len, unsigned cbs, const Sequence& target, Statistics& stats, std::pmr::monotonic_buffer_resource& pool, EMatrixAdjustRule rule) :
    scores(&pool)
{
    TaskTimer timer;
    
    //Masking::get()(target_seq.data(), target_seq.size(), Masking::Algo::SEG);

    const auto c = composition(target);
    scores.resize(32 * AMINO_ACID_COUNT);
    //scores32.resize(32 * AMINO_ACID_COUNT);
    score_min = INT_MAX;
    score_max = INT_MIN;
    array<int, AMINO_ACID_COUNT * AMINO_ACID_COUNT> s;

    //if (cbs == CBS::COMP_BASED_STATS) // || rule == eCompoScaleOldMatrix)
		//throw std::runtime_error("Unsupported CBS code: " + std::to_string(cbs));
        //s = CompositionBasedStats(score_matrix.matrix32_scaled_pointers().data(), query_comp, c, score_matrix.ungapped_lambda(), score_matrix.freq_ratios());
    //else if (cbs == CBS::HAUSER_GLOBAL)
        //throw std::runtime_error("Unsupported CBS code: " + std::to_string(cbs));
        //s = hauser_global(query_comp, c);
    //else
    if (rule == eUserSpecifiedRelEntropy) {
        CompositionMatrixAdjust(query_len, count_true_aa(target), query_comp.data(), c.data(), config.cbs_matrix_scale, score_matrix.ideal_lambda(), score_matrix.joint_probs(), score_matrix.background_freqs(), s, stats);
        stats.inc(Statistics::MATRIX_ADJUST_COUNT, 1);
    }
    else if (rule == eCompoScaleOldMatrix) {
        if (!CompositionBasedStats(score_matrix.matrix32_scaled_pointers().data(), query_comp, c, score_matrix.ungapped_lambda(), score_matrix.freq_ratios(), s)) {
            stats.inc(Statistics::FAILED_COMP_BASED_STATS, 1);
            CompositionMatrixAdjust(query_len, count_true_aa(target), query_comp.data(), c.data(), config.cbs_matrix_scale, score_matrix.ideal_lambda(), score_matrix.joint_probs(), score_matrix.background_freqs(), s, stats);
        }
        else
            stats.inc(Statistics::COMP_BASED_STATS_COUNT, 1);
    }
    else
        throw runtime_error("Unsupported CBS rule: " + std::to_string(rule));

    for (size_t i = 0; i < AMINO_ACID_COUNT; ++i) {
        for (size_t j = 0; j < AMINO_ACID_COUNT; ++j)
            if ((i < 20 || i == MASK_LETTER) && (j < 20 || j == MASK_LETTER)) {
                scores[i * 32 + j] = s[j * AMINO_ACID_COUNT + i];
                if (s[j * AMINO_ACID_COUNT + i] > SCHAR_MAX)
                    s[j * AMINO_ACID_COUNT + i] = SCHAR_MAX;
                else if (s[j * AMINO_ACID_COUNT + i] < SCHAR_MIN)
                    s[j * AMINO_ACID_COUNT + i] = SCHAR_MIN;
                //scores32[i * 32 + j] = s[j * AMINO_ACID_COUNT + i];
                score_min = std::min(score_min, s[j * AMINO_ACID_COUNT + i]);
                score_max = std::max(score_max, s[j * AMINO_ACID_COUNT + i]);
                //std::cerr << s[j * AMINO_ACID_COUNT + i] << ' ';
            }
            else {
                scores[i * 32 + j] = std::max(score_matrix(i, j) * config.cbs_matrix_scale, SCHAR_MIN);
                //scores32[i * 32 + j] = score_matrix(i, j) * config.cbs_matrix_scale;
                //score_min = std::min(score_min, scores32[i * 32 + j]);
                //score_max = std::max(score_max, scores32[i * 32 + j]);
            }
        //std::cerr << std::endl;
    }
    
    stats.inc(Statistics::TIME_MATRIX_ADJUST, timer.microseconds());
    //scores_low = Scores<int8_t>(32, scores.data(), 1, 0, 16);
    //scores_high = Scores<int8_t>(32, scores.data(), 1, 0, 16, 16);
}

}