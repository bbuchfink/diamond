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

#include "chain.h"
#include "../util/util.h"
#include <cstdint>
#include "../util/math/log2_fast.h"

namespace Dna {



void only_keep_best_chains_per_target(std::vector<Chain> &chains, const float chainCutoffPercentage) {

    std::sort(chains.begin(), chains.end(), [](const Chain &a, const Chain &b) {
        return (a.target_id == b.target_id) ? (a.chain_score > b.chain_score) : (a.target_id < b.target_id);
    });

    int last_target_id = chains[0].target_id;
    chains[0].is_primary = 1;
    auto cutoff = int (chainCutoffPercentage * chains[0].chain_score);

    for (size_t i = 1; i < chains.size(); ++i) {
        if (chains[i].chain_score >= cutoff && chains[i].target_id == last_target_id) {
            chains[i].is_primary = 1;
        }
        else if (chains[i].target_id != last_target_id) {
            last_target_id = chains[i].target_id;
            chains[i].is_primary = 1;
            cutoff = int (chainCutoffPercentage * chains[i].chain_score);
        }
    }

    std::sort(chains.begin(), chains.end(), [](const Chain &a, const Chain &b) {
        return a.is_primary > b.is_primary;
    });

    const auto first_secondary = std::find_if(chains.begin(), chains.end(), [](const Chain& chain) {
        return chain.is_primary == 0;
    });

    chains.erase(first_secondary, chains.end());
}


void detect_primary_chains(std::vector<Chain> &chains) {

    std::vector<int> best_secondary_score_of_primary(chains.size(), 0);
    // first chain is always primary
    std::vector<size_t> primary_chain_indices = {0};

    std::vector<int> chain_span;
    chain_span.reserve(chains.size());

    // pre compute chain span
    for (Chain chain : chains) {
        chain_span.push_back(chain.anchors[0].i - chain.anchors.back().i_start());
    }

    for (size_t index_chain = 1; index_chain < chains.size(); ++index_chain) {
        bool isPrimary = true;

        for (size_t index_primary_chain : primary_chain_indices) {
            int overlapLength = chains[index_chain].overlapInQuery(chains[index_primary_chain]);
            if (overlapLength < 1) {
                continue;
            }

            double overlapPercentage = static_cast<double>(overlapLength) / std::min(chain_span[index_chain], chain_span[index_primary_chain]);

            if (overlapPercentage >= MIN_OVERLAP_PERCENTAGE_SECONDARY) {
                isPrimary = false;
                // Chain index_chain is the best secondary to Chain index_primary_chain
                best_secondary_score_of_primary[index_primary_chain] =
                        std::max(best_secondary_score_of_primary[index_primary_chain], chains[index_chain].chain_score);
            }
        }
        if (isPrimary) {
            primary_chain_indices.push_back(index_chain);
        }
    }

    for (size_t i : primary_chain_indices) {
        chains[i].computeMappingQuality(best_secondary_score_of_primary[i]);
    }
}



static int64_t find_chain_start(const int32_t maxDrop, const uint64_t scoreEnd, const uint64_t indexEnd, const AnchorData &anchorData) {

    int64_t index = indexEnd;
    int64_t index_max_score = index;
    int32_t max_score = 0;
    // check for invalid or visited anchor
    if (index < 0 || anchorData.anchor_used[index]) return index;
    // iterate over all anchors in chain while anchor has not been visited
    do {
        // move to the previous anchor in the chain
        index = anchorData.predecessor_anchor[index];
        // compute score difference between anchors (if valid anchor)
        const int32_t score_difference_anchors = index < 0 ? scoreEnd : (int32_t)scoreEnd - anchorData.best_score_anchor[index];
        // new max score? (best extension so far)
        if (score_difference_anchors > max_score) {
            max_score = score_difference_anchors;
            index_max_score = index;
        }
        // chain does not extend if score drops too much between anchors
        else if (max_score - score_difference_anchors > maxDrop) break;
    } while (index > -1 && !anchorData.anchor_used[index]);
    return index_max_score;
}

std::vector<Chain> chain_backtrack(AnchorData &anchorData, const int minChainScore, const int maxDrop, const SeedMatch *seedMatchBegin, bool isReverse) {

    std::vector<ScoreIndexPair<uint64_t, uint64_t>> potential_chain_ends;
    for (int64_t i = 0; i < anchorData.best_score_anchor.size(); ++i)
        if (anchorData.best_score_anchor[i] >= minChainScore)
            potential_chain_ends.emplace_back(anchorData.best_score_anchor[i], i);



    if (potential_chain_ends.empty()) return {};

    if (config.best_hsp_only){

        std::sort(potential_chain_ends.begin(), potential_chain_ends.end(),
        [](const ScoreIndexPair<uint64_t, uint64_t>& a, const ScoreIndexPair<uint64_t, uint64_t>& b) {
            return a.score > b.score; });

        if (potential_chain_ends.size() > 4) {
            potential_chain_ends.resize(4);
        }
    }

    // sort by score
    std::sort(potential_chain_ends.begin(), potential_chain_ends.end());

    std::vector<Chain> detectedChains;
    anchorData.anchor_used.assign(anchorData.anchor_used.size(), false);
    // iterate over all potential end positions (highest to lowest)
    for (int64_t k = potential_chain_ends.size() - 1; k > -1; --k) {

        if (anchorData.anchor_used[potential_chain_ends[k].index]) continue;

        const int64_t chain_start_index = find_chain_start(maxDrop, potential_chain_ends[k].score, potential_chain_ends[k].index, anchorData);


        Chain chain = Chain(isReverse);
        // iterate over all positions in chain, marks used anchors
        int64_t idx_chain_end;
        for (idx_chain_end = potential_chain_ends[k].index; idx_chain_end != chain_start_index; idx_chain_end = anchorData.predecessor_anchor[idx_chain_end]) {
            anchorData.anchor_used[idx_chain_end] = true;
            chain.anchors.emplace_back(seedMatchBegin[idx_chain_end].i(), seedMatchBegin[idx_chain_end].j(), seedMatchBegin[idx_chain_end].ungapped_score());
        }
        // score of current chain
        const int32_t score = idx_chain_end < 0 ? potential_chain_ends[k].score : (int32_t)potential_chain_ends[k].score - anchorData.best_score_anchor[idx_chain_end];
        // valid chain?
        if (score >= minChainScore && !chain.anchors.empty()){
            chain.target_id = seedMatchBegin->id();
            chain.chain_score = score;
            detectedChains.push_back(std::move(chain));
        }
    }
return detectedChains;
}

int32_t compute_score(const SeedMatch &secondMatch, const SeedMatch &firstMatch, const ChainingParameters &chainingParameters) {

    const int32_t distance_query = secondMatch.i_start() - firstMatch.i();
    const int32_t distance_query_end_to_end = secondMatch.i() - firstMatch.i();
    if (distance_query_end_to_end < 1 || distance_query > chainingParameters.MAX_DIST_X)
        return INT32_MIN;

    const int32_t distance_target = secondMatch.j_start() - firstMatch.j();
    const int32_t distance_target_end_to_end = secondMatch.j() - firstMatch.j();
    if (distance_target_end_to_end == 0 || distance_target > chainingParameters.MAX_DIST_Y)
        return INT32_MIN;

    // distance off diagonal (0 is on diagonal)
    const int32_t distance_diagonal = distance_target_end_to_end > distance_query_end_to_end ?
            distance_target_end_to_end - distance_query_end_to_end
            :
            distance_query_end_to_end - distance_target_end_to_end;

    // too big distance on query or target
    if (distance_diagonal > chainingParameters.BAND_WIDTH)
        return INT32_MIN;

    // smaller distance on query or target
    const int32_t distance_skip = std::min(std::abs(distance_target), std::abs(distance_query));
    const int32_t distance_gap_end_to_end = std::min(distance_target_end_to_end, distance_query_end_to_end);

    // initial score: smaller matchSpan or distance
    int32_t score = std::min(secondMatch.ungapped_score(), distance_gap_end_to_end);
    if (distance_diagonal) {
        float lin_pen = chainingParameters.CHAIN_PEN_GAP * (float)distance_diagonal + chainingParameters.CHAIN_PEN_SKIP * (float)distance_skip;
        float log_pen = distance_diagonal > 0 ? log2_approximate(distance_diagonal + 1) : 0.0f;
        score -= (int)(lin_pen + .5f * log_pen);
    }

    return score;
}
  
std::vector<Chain>
chaining_dynamic_program(const ChainingParameters &chainingParameters, const SeedMatch *seedMatchBegin,
                         const SeedMatch *seedMatchEnd, bool isReverse) {
    const auto total_number_matches = std::distance(seedMatchBegin, seedMatchEnd);

    AnchorData anchorData(total_number_matches);
    int32_t total_max_score = 0;

    // fill matrix
    int64_t max_score_index = -1;
    for (int64_t index_second_match = 0; index_second_match < total_number_matches; ++index_second_match) {
        int64_t index_predecessor = -1;
        int32_t max_score = seedMatchBegin[index_second_match].ungapped_score();
        int32_t n_skip = 0;

        int64_t start = 0;
        while (start < index_second_match && seedMatchBegin[index_second_match].j_start() > seedMatchBegin[start].j() + chainingParameters.MAX_DIST_X) ++start;

        start = std::max(start, index_second_match - chainingParameters.MAX_ITERATIONS);

        int64_t index_first_match;
        for (index_first_match = index_second_match - 1; index_first_match >= start; --index_first_match) {
            int32_t score = compute_score(seedMatchBegin[index_second_match], seedMatchBegin[index_first_match],
                                          chainingParameters);
            if (score == INT32_MIN) continue;
            score += anchorData.best_score_anchor[index_first_match];

            if (score > max_score) {
                max_score = score;
                index_predecessor = index_first_match;
                // pending skipped seeds?
                if (n_skip > 0) --n_skip;

            // already in chain?
            } else if (anchorData.pre_predecessor_anchor[index_first_match] == (int32_t)index_second_match) {
                // increment number of skipped seeds and break if too many
                if (++n_skip > chainingParameters.MAX_SKIP) break;
            }
            // updates the seed information if previous seed was part of the chain
            if (anchorData.predecessor_anchor[index_first_match] > -1) anchorData.pre_predecessor_anchor[anchorData.predecessor_anchor[index_first_match]] = index_second_match;
        }
        int64_t end_j = index_first_match;
        if (max_score_index < 0 || seedMatchBegin[index_second_match].j_start() - seedMatchBegin[max_score_index].j() > (int64_t)chainingParameters.MAX_DIST_X) {
            int32_t max = INT32_MIN;
            max_score_index = -1;
            // find seed with the highest score
            for (index_first_match = index_second_match - 1; index_first_match >= start; --index_first_match) {
                if (max < anchorData.best_score_anchor[index_first_match]){
                    max = anchorData.best_score_anchor[index_first_match];
                    max_score_index = index_first_match;
                }
            }
        }
        // valid max score
        if (max_score_index >= 0 && max_score_index < end_j) {
            const int32_t score_extend_max = compute_score(seedMatchBegin[index_second_match],
                                                           seedMatchBegin[max_score_index], chainingParameters);
            // score is valid and higher than the current max score
            if (score_extend_max != INT32_MIN && max_score < score_extend_max + anchorData.best_score_anchor[max_score_index]){
                max_score = score_extend_max + anchorData.best_score_anchor[max_score_index];
                index_predecessor = max_score_index;
            }
        }

        anchorData.best_score_anchor[index_second_match] = max_score;
        anchorData.predecessor_anchor[index_second_match] = index_predecessor;
        anchorData.peak_score_anchor[index_second_match] = index_predecessor > -1 && anchorData.peak_score_anchor[index_predecessor] > max_score ? anchorData.peak_score_anchor[index_predecessor] : max_score;
        if (max_score_index < 0 || (anchorData.best_score_anchor[max_score_index] < anchorData.best_score_anchor[index_second_match]))
            max_score_index = index_second_match;
        // maximum score in entire chaining
        total_max_score = std::max(total_max_score, max_score);
    }

    return chain_backtrack(anchorData, chainingParameters.MIN_CHAIN_SCORE, chainingParameters.BAND_WIDTH, seedMatchBegin, isReverse);

}

}
