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

#include "cbs.h"
#include "../basic/config.h"
#include "score_matrix.h"
#include "../masking/masking.h"

using std::vector;

namespace Stats {

CBS comp_based_stats(0, -1.0, -1.0, -1.0);

CBS::CBS(unsigned code, double query_match_distance_threshold, double length_ratio_threshold, double angle):
    query_match_distance_threshold(-1.0),
    length_ratio_threshold(-1.0),
    angle(50.0)
{
    switch (code) {
    case COMP_BASED_STATS_AND_MATRIX_ADJUST:
        this->angle = 70.0;
        this->query_match_distance_threshold = 0.16;
        this->length_ratio_threshold = 3.0;
    default:
        ;
    }
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
    if (config.comp_based_stats != CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST || a.length() != b.length())
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

TargetMatrix::TargetMatrix(const Composition& query_comp, int query_len, const Sequence& target)
{
    if (!CBS::matrix_adjust(config.comp_based_stats) || target.length() == 0 || query_len == 0)
        return;
    
    vector<Letter> target_seq = target.copy();
    //Masking::get()(target_seq.data(), target_seq.size(), Masking::Algo::SEG);

    //auto c = composition(target);
    auto c = composition(Sequence(target_seq.data(), target_seq.size()));
    EMatrixAdjustRule rule = eUserSpecifiedRelEntropy;
    if (CBS::conditioned(config.comp_based_stats)) {
        rule = s_TestToApplyREAdjustmentConditional(query_len, (int)target.length(), query_comp.data(), c.data(), score_matrix.background_freqs());
        if (rule == eCompoScaleOldMatrix && config.comp_based_stats != CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST)
            return;
    }

    scores.resize(32 * AMINO_ACID_COUNT);
    scores32.resize(32 * AMINO_ACID_COUNT);
    score_min = INT_MAX;
    score_max = INT_MIN;
    vector<int> s;
    
    if (config.comp_based_stats == CBS::COMP_BASED_STATS || rule == eCompoScaleOldMatrix)
        s = CompositionBasedStats(score_matrix.matrix32_scaled_pointers().data(), query_comp, c, score_matrix.ungapped_lambda(), score_matrix.freq_ratios());
    else if (config.comp_based_stats == CBS::HAUSER_GLOBAL)
        s = hauser_global(query_comp, c);
    else
        s = CompositionMatrixAdjust(query_len, count_true_aa(target), query_comp.data(), c.data(), config.cbs_matrix_scale, score_matrix.ideal_lambda(), score_matrix.joint_probs(), score_matrix.background_freqs());
    
    for (size_t i = 0; i < AMINO_ACID_COUNT; ++i) {
        for (size_t j = 0; j < AMINO_ACID_COUNT; ++j)
            if ((i < 20 || i == MASK_LETTER) && (j < 20 || j == MASK_LETTER)) {
                scores[i * 32 + j] = s[j * AMINO_ACID_COUNT + i];
                scores32[i * 32 + j] = s[j * AMINO_ACID_COUNT + i];
                score_min = std::min(score_min, s[j * AMINO_ACID_COUNT + i]);
                score_max = std::max(score_max, s[j * AMINO_ACID_COUNT + i]);
                //std::cerr << s[j * AMINO_ACID_COUNT + i] << ' ';
            }
            else {
                scores[i * 32 + j] = std::max(score_matrix(i, j) * config.cbs_matrix_scale, SCHAR_MIN);
                scores32[i * 32 + j] = score_matrix(i, j) * config.cbs_matrix_scale;
                score_min = std::min(score_min, scores32[i * 32 + j]);
                score_max = std::max(score_max, scores32[i * 32 + j]);
            }
        //std::cerr << std::endl;
    }
}

TargetMatrix::TargetMatrix(const int16_t* query_matrix, const int16_t* target_matrix) :
    scores(32 * AMINO_ACID_COUNT),
    scores32(32 * AMINO_ACID_COUNT),
    score_min(INT_MAX),
    score_max(INT_MIN)
{
    const double f = (double)config.cbs_matrix_scale / (double)CBS::AVG_MATRIX_SCALE;
    for (size_t i = 0; i < AMINO_ACID_COUNT; ++i) {
        for (size_t j = 0; j < AMINO_ACID_COUNT; ++j)
            if (i < 20 && j < 20) {
                const int score = (int)std::round(((double)query_matrix[i * 20 + j] + (double)target_matrix[i * 20 + j]) * f / 2.0);
                scores[i * 32 + j] = score;
                scores32[i * 32 + j] = score;
                score_min = std::min(score_min, score);
                score_max = std::max(score_max, score);
            }
            else {
                const int score = score_matrix(i, j) * config.cbs_matrix_scale;
                scores[i * 32 + j] = score;
                scores32[i * 32 + j] = score;
                score_min = std::min(score_min, score);
                score_max = std::max(score_max, score);
            }
    }
}

}