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

const double MIN_OVERLAP_PERCENTAGE = 0.5;

struct ChainingParameters {
    int MAX_DIST_X = 5000;
    int MAX_DIST_Y = 5000;
    const int BAND_WIDTH = 500;
    const int MAX_SKIP = 25;
    const int MAX_ITER = 5000;
    const int MIN_CHAIN_SCORE = 40;
    const int MIN_NUMBER_MINIMIZERS = 3;
    const float MAP_PERCENTAGE = 0.5;

    const float CHAIN_PEN_GAP;
    const float CHAIN_PEN_SKIP;

    // constructor for additional parameters
    ChainingParameters(float gap, float skip)
            : CHAIN_PEN_GAP(gap), CHAIN_PEN_SKIP(skip) {}
};

struct AnchorData {
    // optimal predecessor for each anchor
    std::vector<int64_t> predecessor_anchor;
    // best score for each anchor
    std::vector<int32_t> best_score_anchor;
    // best score up to now (peak)
    std::vector<int32_t> peak_score_anchor;

    std::vector<int32_t> temp_marking;

    explicit AnchorData(int n) :
            predecessor_anchor(n),
            best_score_anchor(n),
            peak_score_anchor(n),
            temp_marking(n, 0) {}
};



struct Chain {
    explicit Chain(bool rev);
    int32_t chain_score;
    BlockId target_id;
    uint8_t mapping_quality;
    bool reverse;

    struct Anchor {
        Loc i;
        Loc j;
        Anchor(Loc i, Loc j) : i(i), j(j) {}
    };
    // query/target starting position in reverse (i, j)
    std::vector<Anchor> anchors;

    int overlap_query(const Chain &other, const int kmer_size) const {
        return (std::min(anchors[0].i, other.anchors[0].i) + kmer_size) -
                            std::max(anchors.back().i, other.anchors.back().i);

    }

    void compute_mapping_quality(int score_secondary) {
        double sc_ratio = static_cast<double>(score_secondary) / chain_score;
        double quality = 40 * (1 - sc_ratio) * std::min(1.0, static_cast<double>(anchors.size()) / 10) * std::log(chain_score);
        // compress to 0-60
        mapping_quality = static_cast<uint8_t>(quality * 60 / 312);
    }

    bool operator>(const Chain& other)const
    {return this->chain_score > other.chain_score;}

};






std::vector<Chain> chain_dp(int window, int kmer_size, const ChainingParameters &p, const SeedMatch* begin,
                            const SeedMatch* end, bool reverse);
void compute_primary_chains(std::vector<Chain> &chains, int kmer_size);
}

