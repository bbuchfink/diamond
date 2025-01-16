/****
DIAMOND protein aligner
Copyright (C) 2023 Vincent Spath

Code developed by Vincent Spath <vincent.spath@tuebingen.mpg.de>

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

#include "seed_set_dna.h"

namespace Dna {

const double MIN_OVERLAP_PERCENTAGE_SECONDARY = 0.5;

struct ChainingParameters {
    int MAX_DIST_X = 1000; // TODO: should be different for max sensitivity and assembly alignment (up to 5000 for less kmers)
    int MAX_DIST_Y = 1000;
    const int BAND_WIDTH = 300; // TODO: should be different for max sensitivity and assembly alignment (500)
    const int MAX_SKIP = 25;
    const int MAX_ITERATIONS = 3000; // TODO: adapt all for possible short mode? band usw abhängig machen von der query length?
    const float MAP_PERCENTAGE_TARGET = 0.99f;

    const int MIN_CHAIN_SCORE;
    const float CHAIN_PEN_GAP;
    const float CHAIN_PEN_SKIP;
    const float MAX_OVERLAP_EXTENSION;

    ChainingParameters(float gap, float skip, int min_chain_score, float max_overlap_extension)
            : MIN_CHAIN_SCORE(min_chain_score), CHAIN_PEN_GAP(gap), CHAIN_PEN_SKIP(skip), MAX_OVERLAP_EXTENSION(max_overlap_extension) {}
};

struct AnchorData {
    std::vector<int64_t> predecessor_anchor;
    std::vector<int32_t> best_score_anchor;
    std::vector<int32_t> peak_score_anchor;
    std::vector<int32_t> pre_predecessor_anchor;
    std::vector<bool> anchor_used;

    explicit AnchorData(size_t n) :
            predecessor_anchor(n),
            best_score_anchor(n),
            peak_score_anchor(n),
            pre_predecessor_anchor(n, 0),
            anchor_used(n){}
};



struct Chain {
    explicit Chain(bool rev);
    int32_t chain_score;
    BlockId target_id;
    uint8_t mapping_quality;
    uint8_t is_primary = 0;
    bool reverse;

    struct Anchor {
        Loc i;
        Loc j;
        int span;
        Loc i_start() const {
            return i - span;
        }
        Loc j_start() const {
            return j - span;
        }
        Anchor(Loc i, Loc j, int span) : i(i), j(j), span(span) {}
    };
    // query/target starting position in reverse (i, j)
    std::vector<Anchor> anchors;

    explicit Chain(bool rev) : chain_score(0), target_id(0), mapping_quality(0), reverse(rev), anchors() {}

    int overlapInQuery(const Chain &otherChain) const {
        return (std::min(anchors[0].i, otherChain.anchors[0].i)) -
               std::max(anchors.back().i_start(), otherChain.anchors.back().i_start());
    }

    int overlapInTarget(const Chain &otherChain) const {
        return (std::min(anchors[0].j, otherChain.anchors[0].j)) -
               std::max(anchors.back().j_start(), otherChain.anchors.back().j_start());
    }

    void computeMappingQuality(int scoreSecondaryChain) {
        double score_ratio = static_cast<double>(scoreSecondaryChain) / chain_score;
        double quality_score = 40 * (1 - score_ratio) * std::min(1.0, static_cast<double>(anchors.size()) / 10) * std::log(chain_score);
        // compress to 0-60
        mapping_quality = static_cast<uint8_t>(quality_score * 60 / 312);
    }

    bool operator>(const Chain& otherChain)const
    {return this->chain_score > otherChain.chain_score;}
};


template <typename T1, typename T2>
struct ScoreIndexPair {
    static_assert(std::is_arithmetic<T1>::value, "T1 must be an arithmetic type");
    static_assert(std::is_arithmetic<T2>::value, "T2 must be an arithmetic type");
    T1 score;
    T2 index;

    ScoreIndexPair() : score(T1()), index(T2()) {}
    ScoreIndexPair(const T1& first, const T2& second) : score(first), index(second) {}

    // sort ascending
    bool operator<(const ScoreIndexPair& other) const {
        if (score != other.score)
            return score < other.score;
        else
            return index < other.index;
    }
};




std::vector<Chain>
chaining_dynamic_program(const ChainingParameters &chainingParameters, const SeedMatch *seedMatchBegin,
                         const SeedMatch *seedMatchEnd, bool isReverse);
void detect_primary_chains(std::vector<Chain> &chains);
void only_keep_best_chains_per_target(std::vector<Chain> &chains, float chainCutoffPercentage);

}