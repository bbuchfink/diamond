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

#include "extension_chain.h"
#include "seed_set_dna.h"
#include "timer.h"
#include "chain.h"
#include "alignment.h"
#include "extension_seed_matches.h"

namespace Dna {

int compute_residue_matches_of_chain(const std::vector<Chain::Anchor> &anchors, const int kmer_size) {

    int totalMatchingResidues = anchors.back().span;

    for (size_t i = anchors.size() - 1; i > 0; --i) {
        totalMatchingResidues += std::min(std::min(anchors[i-1].span, anchors[i-1].i - anchors[i].i),
                                          anchors[i-1].j - anchors[i].j);
    }

    return totalMatchingResidues;
}

Extension::Match build_map_hsp(const Search::Config &cfg, const BlockId targetBlockId, const Chain *beginChain, const Chain *endChain) {

    Extension::Match match = Extension::Match(targetBlockId, cfg.target->seqs()[targetBlockId], ::Stats::TargetMatrix(), 0, 0);

    for (auto chain = beginChain; chain != endChain; ++chain) {
        //if (!chain->is_primary) continue; // only map primary chains
        Hsp mapHsp = Hsp();

        mapHsp.query_range.begin_ = chain->anchors.back().i_start();
        mapHsp.subject_range.begin_ = chain->anchors.back().j_start();
        mapHsp.query_range.end_ = chain->anchors[0].i;
        mapHsp.subject_range.end_ = chain->anchors[0].j;

        mapHsp.identities = compute_residue_matches_of_chain(chain->anchors, shapes[0].length_);
        mapHsp.length = std::max(chain->anchors[0].i - chain->anchors.back().i_start(),
                              chain->anchors[0].j - chain->anchors.back().j_start()); // approx, real alignment length might be bigger
        mapHsp.mapping_quality = chain->mapping_quality;
        mapHsp.n_anchors = chain->anchors.size();

        mapHsp.transcript.push_terminator();
        mapHsp.target_seq = cfg.target->seqs()[targetBlockId];
        mapHsp.query_source_range = mapHsp.query_range;
        mapHsp.subject_source_range = chain->reverse ? Interval(mapHsp.subject_range.end_, mapHsp.subject_range.begin_) : Interval(
                mapHsp.subject_range.begin_, mapHsp.subject_range.end_);
        mapHsp.frame = chain->reverse + 2;

        match.hsp.push_back(std::move(mapHsp));
    }

    return match;
}

std::pair<Hsp, int>
extend_new_at_peak(const Search::Config &cfg, Cigar &extension, const Chain *chain, const Sequence &query,
                   const Sequence &target, const BlockId targetBlockId, const int start_i, const int start_j,
                   const int now_anchor_index) {
    auto erase_it = extension.getCigarDataConst().begin() + extension.peakScoreCigarIndex;
    extension.getCigarData().erase(erase_it, extension.getCigarDataConst().end());
    extension.score = extension.peakScore;

    // extend only max to next anchor? TODO: delete commented code when sure
    /*Sequence query_right = query.subseq(chain->anchors[extension.peakScoreAnchorIndex].i,
                               chain->anchors[extension.peakScoreAnchorIndex - 1].i_start());
    Sequence target_right = target.subseq(chain->anchors[extension.peakScoreAnchorIndex].j,
                                 chain->anchors[extension.peakScoreAnchorIndex - 1].j_start());
*/
    // if query or target overlap just align until end
    //if (query_right.length() < 1 || target_right.length() < 1) {
    const auto j_end = chain->anchors[extension.peakScoreAnchorIndex].j;
    Sequence query_right = query.subseq(chain->anchors[extension.peakScoreAnchorIndex].i, query.length());
    Sequence target_right = target.subseq(j_end,
                                     std::min(target.length(),
                                                       std::min(j_end + config.band_extension + (int) query_right.length(),
                                                                j_end + (int) query_right.length() * 2)));
    //}

    config.dna_extension == DNAExtensionAlgo::WFA ?
    compute_wfa_cigar(cfg, query_right.to_string(), extension, false, false, target_right.to_string(), WFA_BAND_EXTENSION)
                                                  :
    compute_ksw_cigar(target_right, query_right, cfg, KSW_FLAG_R, extension, config.zdrop_extension, config.band_extension);

    Hsp hsp = build_hsp_from_cigar(extension, cfg.target->seqs()[targetBlockId], query, start_i,
                                   start_j,
                                   chain->reverse, cfg);

    // TODO: return which anchor (as heuristic) for next extension? After peak or from now on? Save multiple peak scores?
    //return std::make_pair(hsp, extension.peakScoreAnchorIndex - 1);
    return std::make_pair(hsp, now_anchor_index);
}


std::pair<Hsp,int> extend_between_anchors(const Search::Config &cfg, const BlockId targetBlockId, const Chain *chain, const Sequence &query, const Sequence &target, int anchorIdx) {

    // start position of chain
    const int start_i = chain->anchors[anchorIdx].i_start();
    const int start_j = chain->anchors[anchorIdx].j_start();

    Cigar extension(anchorIdx * 3);

    if (start_i > 0 && start_j > 0) {
        int prev_anchor_i = 0;
        int prev_anchor_j = 0;

        // extend maximal to end anchor of last alignment
        if (anchorIdx != static_cast<int>(chain->anchors.size()) - 1) {
            if (chain->anchors[anchorIdx + 1].i <= start_i && chain->anchors[anchorIdx + 1].j <= start_j) {
                prev_anchor_i = chain->anchors[anchorIdx + 1].i;
                prev_anchor_j = chain->anchors[anchorIdx + 1].j;
            }
        }


        // extend to the left at start (until previous anchor if dropped)
        std::vector<Letter> query_left = query.subseq(prev_anchor_i, start_i).reverse();
        std::vector<Letter> target_left = target.subseq(std::max(prev_anchor_j,
                                                                std::max(start_j - (int)query_left.size() - config.band_extension,
                                                                        start_j - ((int)query_left.size() * 2))),
                                                        start_j).reverse();

            config.dna_extension == DNAExtensionAlgo::WFA ?
            compute_wfa_cigar(cfg, ((Sequence) query_left).to_string(), extension, true, false,
                            ((Sequence) target_left).to_string(), WFA_BAND_EXTENSION)
                                                                        :
            compute_ksw_cigar(target_left, query_left, cfg, KSW_FLAG_L, extension, config.zdrop_extension, config.band_extension);
    }
    else {
        extension.setMaxValues(-1, -1);
    }

    int anchor_distance_query= INT_MAX;
    int anchor_distance_target = INT_MAX;
    // iterate anchors from start to end (reverse ordered)
    for (; anchorIdx > 0; --anchorIdx) {

        const auto span_current_anchor = chain->anchors[anchorIdx].span;
        const auto span_next_anchor = chain->anchors[anchorIdx - 1].span;

        // anchor no overlap with previous anchor
        if (anchor_distance_query > span_current_anchor && anchor_distance_target > span_current_anchor) {
            extension.extendCigar(span_current_anchor, 'M');
            extension.score += span_current_anchor * cfg.score_builder->reward();
        }

        if (extension.score > extension.peakScore) {
            extension.peakScore = extension.score;
            extension.peakScoreCigarIndex = extension.getCigarDataConst().size();
            extension.peakScoreAnchorIndex = anchorIdx;
        }


        anchor_distance_query = chain->anchors[anchorIdx - 1].i - chain->anchors[anchorIdx].i;
        anchor_distance_target = chain->anchors[anchorIdx - 1].j - chain->anchors[anchorIdx].j;

        // Case 1: Anchors don't overlap in Query and Target
        if (anchor_distance_query > span_next_anchor && anchor_distance_target > span_next_anchor) {
            // extend to the right between anchors
            Sequence query_right = query.subseq(chain->anchors[anchorIdx].i, chain->anchors[anchorIdx - 1].i_start());
            Sequence target_right = target.subseq(chain->anchors[anchorIdx].j, chain->anchors[anchorIdx - 1].j_start());

            const int alignment_band = std::abs(query_right.length() - target_right.length())
                    + std::min(config.band_global, (std::min(query_right.length(), target_right.length())) / 2);

            auto z_dropped = (config.dna_extension == DNAExtensionAlgo::WFA ?
                    compute_wfa_cigar(cfg, query_right.to_string(), extension, false, true,
                                                target_right.to_string(), alignment_band)
                                                :
                    compute_ksw_cigar(target_right, query_right, cfg, KSW_FLAG_G, extension,
                                      config.zdrop_global, alignment_band));

            // finish alignment and start new one
            if (z_dropped == DROPPED || z_dropped == NEGATIVE_SCORE) {

                if (extension.score >= extension.peakScore) {
                    Hsp hsp = build_hsp_from_cigar(extension, cfg.target->seqs()[targetBlockId], query, start_i,
                                                   start_j,
                                                   chain->reverse, cfg);
                    return std::make_pair(hsp, anchorIdx - 1);
                }
                else {
                    return extend_new_at_peak(cfg, extension, chain, query, target, targetBlockId, start_i, start_j, anchorIdx - 1);
                }
            }
        }


        // Case 2: Anchors overlap more in Query than in Target
        else if (anchor_distance_query <= span_next_anchor && anchor_distance_query < anchor_distance_target) {
            // push deletions
            const int number_of_gaps = anchor_distance_target - anchor_distance_query;
            extension.extendCigar(number_of_gaps, 'D');
            extension.score -= (number_of_gaps * cfg.score_builder->gap_extend()) + cfg.score_builder->gap_open();

            // push matches for last part of second anchor
            extension.extendCigar(anchor_distance_query, 'M');
            extension.score += anchor_distance_query * cfg.score_builder->reward();
        }

        // Case 3: Anchors overlap more in Target than in Query
        else if (anchor_distance_target <= span_next_anchor && anchor_distance_target < anchor_distance_query) {
            // push insertions
            const int number_of_gaps = anchor_distance_query - anchor_distance_target;
            extension.extendCigar(number_of_gaps, 'I');
            extension.score -= (number_of_gaps * cfg.score_builder->gap_extend()) + cfg.score_builder->gap_open();

            // push matches for last part of second anchor
            extension.extendCigar(anchor_distance_target, 'M');
            extension.score += anchor_distance_target * cfg.score_builder->reward();
        }
        else {
            throw std::runtime_error("Error in chaining extension: No anchor overlap case matched");
        }
    }

    const auto span_last_anchor = chain->anchors.front().span;
    // no overlap with last anchor
    if (anchor_distance_query > span_last_anchor && anchor_distance_target > span_last_anchor) {
        extension.extendCigar(span_last_anchor, 'M');
        extension.score += span_last_anchor * cfg.score_builder->reward();
    }

    const auto j_end_last_anchor = chain->anchors.front().j;
    Sequence query_right = query.subseq(chain->anchors.front().i, query.length());
    Sequence target_right = target.subseq(j_end_last_anchor,
                                          std::min(target.length(),
                                                   std::min(j_end_last_anchor + config.band_extension + (int) query_right.length(),
                                                            j_end_last_anchor + (int) query_right.length() * 2)));

    config.dna_extension == DNAExtensionAlgo::WFA ?
    compute_wfa_cigar(cfg, query_right.to_string(), extension, false, false, target_right.to_string(), WFA_BAND_EXTENSION)
    :
    compute_ksw_cigar(target_right, query_right, cfg, KSW_FLAG_R, extension, config.zdrop_extension, config.band_extension);

    if (extension.score >= extension.peakScore) {
        Hsp hsp = build_hsp_from_cigar(extension, cfg.target->seqs()[targetBlockId], query, start_i, start_j,
                                       chain->reverse, cfg);

        return std::make_pair(hsp, anchorIdx - 1);
    }
    else {
        return extend_new_at_peak(cfg, extension, chain, query, target, targetBlockId, start_i, start_j, anchorIdx - 1);
    }
}



Extension::Match
extend_chains(const Search::Config &cfg, const BlockId targetBlockId, const Chain *beginChain, const Chain *endChain,
              const Sequence &query, const Sequence &query_reverse, ChainingParameters &chainingParameters) {

    Extension::Match match = Extension::Match(targetBlockId, cfg.target->seqs()[targetBlockId], ::Stats::TargetMatrix(), 0, 0);

    const Sequence &target = cfg.target->seqs()[targetBlockId];


    ExtensionTimer timer;
    auto start_extend = std::chrono::high_resolution_clock::now();

    std::vector<std::tuple<int, int, int, int>> range_extended_chains;



    for (auto chain = beginChain; chain != endChain; ++chain) {

        // do not extend chains that overlap with percentage previous extensions
        const int query_range = chain->anchors[0].i - chain->anchors.back().i_start();
        const int target_range = chain->anchors[0].j - chain->anchors.back().j_start();
        bool do_not_extend = false;

        for (const auto& [extQueryStart, extQueryEnd, extTargetStart, extTargetEnd] : range_extended_chains) {
            const int query_overlap = std::min(chain->anchors[0].i, extQueryEnd) - std::max(chain->anchors.back().i_start(), extQueryStart);
            const int target_overlap = std::min(chain->anchors[0].j, extTargetEnd) - std::max(chain->anchors.back().j_start(), extTargetStart);

            // TODO: vielleicht im Verhältnis zur Länge der kleineren span?  -> bereits alignierte chains hatten größeren Score, falls abgebrochen, dann trotzdem sonnvoll
            if (query_overlap > int(chainingParameters.MAX_OVERLAP_EXTENSION * query_range)
            && target_overlap > int(chainingParameters.MAX_OVERLAP_EXTENSION * target_range)) {
                do_not_extend = true;
                break;
            }
        }
        if (do_not_extend) {
            continue;
        }


        int anchor_idx = static_cast<int>(chain->anchors.size()) - 1;
        do {
            auto [hsp_result, current_index] = chain->reverse ?
                    extend_between_anchors(cfg, targetBlockId, chain, query_reverse, target, anchor_idx)
                    :
                    extend_between_anchors(cfg, targetBlockId, chain, query, target, anchor_idx);


            if (hsp_result.evalue < config.max_evalue) {
                // TODO: put this maybe in front of the condition?
                range_extended_chains.emplace_back(hsp_result.query_range.begin_, hsp_result.query_range.end_,
                                                   hsp_result.subject_range.begin_, hsp_result.subject_range.end_);
                match.hsp.push_back(std::move(hsp_result));
            }

            anchor_idx = current_index;
        } while (anchor_idx > -1);

    }


    auto end_extend = std::chrono::high_resolution_clock::now();
    timer.update(4,end_extend-start_extend);
        {
            std::lock_guard<std::mutex> lock(cfg.timer->mtx);
            *(cfg.timer) += timer;
        }

    if(config.best_hsp_only) {
        auto best_hsp_it = std::max_element(match.hsp.begin(), match.hsp.end(), [](const Hsp &a, const Hsp &b) {
            return a.score < b.score;
        });

        if (best_hsp_it != match.hsp.end()) {
            Hsp best_hsp = std::move(*best_hsp_it);
            match.hsp.clear();
            match.hsp.push_back(std::move(best_hsp));
        }
    }
    

    return match;
}



std::vector<Chain> compute_chains(const Search::Config &cfg, const Sequence &query, bool isReverse, const ChainingParameters &chainingParameters) {

        auto time_start_seed = std::chrono::high_resolution_clock::now();

    auto seed_hits = seed_lookup(query, cfg.target->seqs(), cfg.dna_ref_index.get(), cfg.minimizer_window);

        auto time_end_seed = std::chrono::high_resolution_clock::now();

        ExtensionTimer timer_seed;
        timer_seed.update(5, time_end_seed - time_start_seed);
        {
        std::lock_guard<std::mutex> lock(cfg.timer->mtx);
        *(cfg.timer) += timer_seed;
        }

    if (seed_hits.empty()) {
        return {};
    }


        auto time_start_ungapped = std::chrono::high_resolution_clock::now();


    seed_hits = merge_and_extend_seeds(seed_hits, query, cfg);


        auto time_end_ungapped = std::chrono::high_resolution_clock::now();

        ExtensionTimer timer_ungapped;
        timer_seed.update(1, time_end_ungapped - time_start_ungapped);
        {
            std::lock_guard<std::mutex> lock(cfg.timer->mtx);
            *(cfg.timer) += timer_seed;
        }


    // Sort SeedMatch vector by BlockId and target location j
    std::sort(seed_hits.begin(), seed_hits.end(), [](const SeedMatch &a, const SeedMatch &b) {
        return a.id() < b.id() || (a.id() == b.id() && a.j() < b.j());
    });

    auto hits_it = merge_keys(seed_hits.begin(), seed_hits.end(), [](const SeedMatch &hit1) { return hit1.id(); });

        auto time_start_chain = std::chrono::high_resolution_clock::now();


    std::vector<Chain> chains;
    // process subsets based on unique target id
    while (hits_it.good()) {

        auto new_chains = chaining_dynamic_program(chainingParameters, hits_it.begin().operator->(),
                                                   hits_it.end().operator->(), isReverse);

        chains.insert(chains.end(), std::make_move_iterator(new_chains.begin()),
                      std::make_move_iterator(new_chains.end()));

        ++hits_it;
    }

        auto time_end_chain = std::chrono::high_resolution_clock::now();

        ExtensionTimer timer_chain;
        timer_chain.update(6, time_end_chain - time_start_chain);

        {
        std::lock_guard<std::mutex> lock(cfg.timer->mtx);
        *(cfg.timer) += timer_chain;
        }

    return chains;

}

std::vector<Extension::Match> chaining_and_extension(const Search::Config &cfg, const Sequence &query, const Sequence &query_reverse) {

    std::vector<Extension::Match> matches;

    ChainingParameters chainingParameters(static_cast<float>(cfg.chain_pen_gap), static_cast<float>(cfg.chain_pen_skip), cfg.min_chain_score, static_cast<float>(cfg.max_overlap_extension)); // static?

    auto chains = compute_chains(cfg, query, false, chainingParameters);
    auto chains_reverse = compute_chains(cfg, query_reverse, true, chainingParameters);

    if (chains.empty() && chains_reverse.empty()) {
        return matches;
    }

    chains.insert( chains.end(), std::make_move_iterator(chains_reverse.begin()),
                   std::make_move_iterator(chains_reverse.end()));


    //if (config.best_hsp_only) {
    //    only_keep_best_chains_per_target(chains, chainingParameters.MAP_PERCENTAGE_TARGET);
    //}

    // sort chains by score
    std::sort(chains.begin(), chains.end(), std::greater<>());

    // compute chain quality
    if (config.chaining_out)
        detect_primary_chains(chains);

    // filter chains by score threshold // TODO: maybe remove from sensitivity modes and only filter when threshold was set (!= 0)
    const int map_score_threshold = chains[0].chain_score * cfg.chain_fraction_align;
    auto lower_bound = std::lower_bound(chains.begin(), chains.end(), map_score_threshold, [](const Chain& chain, int threshold) {
        return chain.chain_score >= threshold;
    });
    chains.erase(lower_bound, chains.end());

    std::sort(chains.begin(), chains.end(), [](const Chain &a, const Chain &b) {
        return (a.target_id == b.target_id) ? (a.chain_score > b.chain_score) : (a.target_id < b.target_id);
    });

    auto iterator_chains_target = merge_keys(chains.begin(), chains.end(),
                                             [](const Chain &chain1) { return chain1.target_id; });

    while (iterator_chains_target.good()) {
        Extension::Match match = config.chaining_out ?
                                 build_map_hsp(cfg, iterator_chains_target.key(), iterator_chains_target.begin().operator->(),
                                               iterator_chains_target.end().operator->()) :
                                 extend_chains(cfg, iterator_chains_target.key(),
                                               iterator_chains_target.begin().operator->(),
                                               iterator_chains_target.end().operator->(), query, query_reverse,
                                               chainingParameters);

        if (!match.hsp.empty()) {
            matches.push_back(std::move(match));
        }

        ++iterator_chains_target;
    }

    return matches;
}

}
