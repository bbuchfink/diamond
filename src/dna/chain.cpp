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


 /**
  * find the start of a chain
  * @param max_drop
  * @param score_end
  * @param index_end
  * @param aD
  * @return index of the start of the chain
  */
static int64_t chain_start(const int32_t max_drop, const uint64_t score_end, const uint64_t index_end, const AnchorData &aD) {

    int64_t i = index_end;
    int64_t max_i = i;
    int32_t max_s = 0;
    // check for invalid or visited anchor
    if (i < 0 || aD.temp_marking[i] != 0) return i;
    // iterate over all anchors in chain while anchor has not been visited
    do {
        // move to the previous anchor in the chain
        i = aD.predecessor_anchor[i];
        // compute score difference between anchors (if valid anchor)
        const int32_t s = i < 0 ? score_end : (int32_t)score_end - aD.best_score_anchor[i];
        // new max score? (best extension so far)
        if (s > max_s) {
            max_s = s;
            max_i = i;
        }
        // chain does not extend if score drops too much between anchors
        else if (max_s - s > max_drop) break;
    } while (i >= 0 && aD.temp_marking[i] == 0);
    return max_i;
}


/**
 * backtrack to find chains
 * @param n
 * @param aD
 * @param min_chain_score
 * @param max_drop
 * @param hits
 * @return vector of chains
 */
std::vector<Chain> chain_backtrack(AnchorData &aD, const int min_chain_score, const int max_drop, const SeedMatch *begin, bool reverse) {

    // potential starting positions of valid chains (score, index)
    std::vector<std::pair<uint64_t, uint64_t>> z_end;
    for (int64_t i = 0; i < aD.best_score_anchor.size(); ++i)
        if (aD.best_score_anchor[i] >= min_chain_score) z_end.emplace_back(aD.best_score_anchor[i], i);

    // sort by score
    std::sort(z_end.begin(), z_end.end());

    // size of z_end is the number of potential end positions of valid chains
    const int64_t n_z = z_end.size();
    if (n_z == 0) return {};

    // chains found during the chaining process
    std::vector<Chain> chains;
    aD.temp_marking.assign(aD.temp_marking.size(), 0);
    // iterate over all potential end positions (highest to lowest)
        for (int64_t k = n_z - 1; k >= 0; --k) {
            // position not been used yet?
            if (aD.temp_marking[z_end[k].second] != 0) continue;
            // calculate the end of the current chain
            const int64_t start_i = chain_start(max_drop, z_end[k].first, z_end[k].second, aD);
            // iterate over all positions in chain, marks used anchors
            Chain chain_t = Chain(reverse);
            int64_t i; // so? brauche es ja für sc später
            for (i = z_end[k].second; i != start_i; i = aD.predecessor_anchor[i]) {
                aD.temp_marking[i] = 1;
                chain_t.anchors.emplace_back(begin[i].i(), begin[i].j());
            }
            // score of current chain
            const int32_t sc = i < 0 ? z_end[k].first : (int32_t)z_end[k].first - aD.best_score_anchor[i];
            // valid chain?
            // prev_n_chain_anchors stores the number of anchors before processing the current chain
            if (sc >= min_chain_score && !chain_t.anchors.empty()){
                chain_t.target_id = begin->id();
                chain_t.chain_score = sc;
                chains.push_back(std::move(chain_t));
            }
        }
    return chains;
}


/**
 * compute score between two anchors
 * @param hit_i
 * @param hit_j
 * @param q_span
 * @param p
 * @return score of the extension
 */
int32_t compute_score(const SeedMatch &hit_i, const SeedMatch &hit_j, int q_span, const ChainingParameters &p) {

    // distance on query
    const int32_t dq = hit_i.i() - hit_j.i();
    if (dq <= 0 || dq > p.MAX_DIST_X) return INT32_MIN;

    // distance on target
    const int32_t dr = hit_i.j() - hit_j.j();
    if (dr == 0 || dr > p.MAX_DIST_Y) return INT32_MIN;

    // absolute difference (in positions) between query and target
    const int32_t dd = dr > dq? dr - dq : dq - dr;

    // too big distance on query or target
    if (dd > p.BAND_WIDTH) return INT32_MIN;

    // smaller distance on query or target (gap)
    const int32_t dg = std::min(dr, dq);

    // initial score: smaller q_span or gap
    int32_t sc = std::min(q_span, dg);
    if (dd || dg > q_span) {
        float lin_pen = p.CHAIN_PEN_GAP * (float)dd + p.CHAIN_PEN_SKIP * (float)dg;
        float log_pen = dd >= 1? log2_ap(dd + 1) : 0.0f; // log2() only works for dd>=2
        sc -= (int)(lin_pen + .5f * log_pen);
    }

    return sc;
}



/**
* identifies the primary chains and their mapping quality
* @param chains
* @param kmer_size
*/
 void compute_primary_chains(std::vector<Chain> &chains, const int kmer_size) { // new version

     std::vector<int> score_secondary(chains.size(), 0);
     // first chain is always primary
     std::vector<size_t> primary_chains = {0};
     std::vector<int> chain_span;
     chain_span.reserve(chains.size());
     // pre compute chain span
     for (Chain chain : chains) {
         chain_span.push_back(chain.anchors[0].i + kmer_size - chain.anchors.back().i);
     }

     // iterate over all chains
     for (size_t i = 1; i < chains.size(); ++i) {
         // chain overlaps with a primary?
         bool primary = true;
         for (size_t c : primary_chains) {
             // calculate overlap length between chains
             int overlapLength = chains[i].overlap_query(chains[c], kmer_size);
             // no overlap
             if (overlapLength <= 0) {
                 continue;
             }

             // calculate overlap percentage of shorter chain
             double overlapPercentage = static_cast<double>(overlapLength) / std::min(chain_span[i], chain_span[c]);

             if (overlapPercentage >= MIN_OVERLAP_PERCENTAGE) {
                 primary = false;
                 // Chain i is the best secondary to Chain c
                 score_secondary[c] = std::max(score_secondary[c], chains[i].chain_score);
             }
         }
         // primary chain added to vector
         if (primary) {
             primary_chains.push_back(i);
         }
     }

     // mapping quality
     for (size_t i : primary_chains) {
         chains[i].compute_mapping_quality(score_secondary[i]);
     }
}


/**
 * dynamic programming chaining algorithm
 * @param hits
 * @param window
 * @param kmer_size
 * @param p
 * @return vector of chains
 */
std::vector<Chain> chain_dp(const int window, const int kmer_size, const ChainingParameters &p, const SeedMatch* begin,
                            const SeedMatch* end, bool reverse) {

    // number of total matches
    const auto n = std::distance(begin, end);

    const int32_t max_drop = p.BAND_WIDTH; // zur Übersicht und testen
    int32_t mmax_f = 0;

    // initialize vectors for chaining
    AnchorData aD(n);

    //code is for 1 vector of matches for 1 query and 1 target
    int64_t max_i_s = -1;
    for (int64_t i = 0; i < n; ++i) {
        int64_t max_j = -1;
        int32_t max_f = kmer_size;
        int32_t n_skip = 0;
        // increase st until same id or t_pos of st is in range of i
        int64_t st = 0;
        while (st < i &&  begin[i].j() > begin[st].j() + p.MAX_DIST_X) ++st;
        // stay in range of max iterations
        st = std::max(st, i - p.MAX_ITER);
        // iterate over all hits from st to i
        int64_t j;
        for (j = i - 1; j >= st; --j) {
            int32_t sc = compute_score(begin[i], begin[j], kmer_size, p);
            if (sc == INT32_MIN) continue;
            sc += aD.best_score_anchor[j];
            // new max score?
            if (sc > max_f) {
                max_f = sc;
                max_j = j;
                // pending skipped seeds?
                if (n_skip > 0) --n_skip;
            // already in chain?
            } else if (aD.temp_marking[j] == (int32_t)i) {
                // increment number of skipped seeds and break if too many
                if (++n_skip > p.MAX_SKIP) break;
            }
            // updates the seed information if previous seed was part of the chain
            if (aD.predecessor_anchor[j] >= 0) aD.temp_marking[aD.predecessor_anchor[j]] = i;
        }
        int64_t end_j = j;
        if (max_i_s < 0 || begin[i].j() - begin[max_i_s].j() > (int64_t)p.MAX_DIST_X) {
            int32_t max = INT32_MIN;
            max_i_s = -1;
            // find seed with the highest score
            for (j = i - 1; j >= st; --j) {
                if (max < aD.best_score_anchor[j]){
                    max = aD.best_score_anchor[j];
                    max_i_s = j;
                }
            }
        }
        // valid max score
        if (max_i_s >= 0 && max_i_s < end_j) {
            // score of extending the current anchor to the best scoring anchor
            const int32_t tmp = compute_score(begin[i], begin[max_i_s], kmer_size, p);
            // score is valid and higher than the current max score
            if (tmp != INT32_MIN && max_f < tmp + aD.best_score_anchor[max_i_s]){
                max_f = tmp + aD.best_score_anchor[max_i_s];
                max_j = max_i_s;
            }
        }
        // setting max score at seed i and index of the seed that contributes to the maximum score of the chain ending at seed i
        // peak_score_anchor keeps the peak score up to i; best_scores_anchors is the score ending at i, not always the peak
        aD.best_score_anchor[i] = max_f;
        aD.predecessor_anchor[i] = max_j;
        aD.peak_score_anchor[i] = max_j >= 0 && aD.peak_score_anchor[max_j] > max_f? aD.peak_score_anchor[max_j] : max_f;
        if (max_i_s < 0 || (aD.best_score_anchor[max_i_s] < aD.best_score_anchor[i]))
            max_i_s = i;
        // maximum score in entire chaining
        mmax_f = std::max(mmax_f, max_f);
    }

    return chain_backtrack(aD, p.MIN_CHAIN_SCORE, max_drop, begin, reverse);

}

// Chain constructor
    Chain::Chain(bool rev) : chain_score(0), anchors(), target_id(0), mapping_quality(0), reverse(rev) {}

}
