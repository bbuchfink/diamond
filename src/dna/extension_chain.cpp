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
#include "extension.h"
#include "seed_set_dna.h"
#include "timer.h"
#include "chain.h"



namespace Dna {


    /**
 * Computes the number of matching residues of all anchors in a chain
 * @param anchors
 * @param kmer_size
 * @return number of matching residues
 */
int compute_residue_matches(const std::vector<Chain::Anchor> &anchors, const int kmer_size) {
    // last one is always kmer_size
    int totalMatchingResidues = kmer_size;

    for (size_t i = anchors.size() - 1; i > 0; --i) {
        totalMatchingResidues += std::min(std::min(kmer_size, anchors[i-1].i- anchors[i].i),
                                          anchors[i-1].j - anchors[i].j);
    }

    return totalMatchingResidues;
}


/**
 * Builds a match object of HSPs from all chains of a target
 * @param cfg
 * @param id
 * @param query
 * @param begin
 * @param end
 * @param reverse
 * @return match object
 */
Extension::Match build_map_HSP(const Search::Config &cfg, const BlockId id, const Chain *begin, const Chain *end) {
    Extension::Match m = Extension::Match(id, cfg.target->seqs()[id], ::Stats::TargetMatrix(), 0, 0);
    for (auto chain = begin; chain != end; ++chain) {
        //if (!chain->is_primary) continue; // only map primary chains
        Hsp out = Hsp();

        out.query_range.begin_ = chain->anchors.back().i;
        out.subject_range.begin_ = chain->anchors.back().j;
        out.query_range.end_ = chain->anchors[0].i + shapes[0].length_; // minimap scheint nicht 1 abzuziehen in Ausgabe
        out.subject_range.end_ = chain->anchors[0].j + shapes[0].length_;

        out.identities = compute_residue_matches(chain->anchors, shapes[0].length_);
        out.length = std::max(chain->anchors[0].i - chain->anchors.back().i + shapes[0].length_,
                              chain->anchors[0].j - chain->anchors.back().j + shapes[0].length_); // Block size: better approximation possible
        out.mapping_quality = chain->mapping_quality;
        out.n_anchors = chain->anchors.size();

        out.transcript.push_terminator();
        out.target_seq = cfg.target->seqs()[id];
        out.query_source_range = out.query_range;
        out.subject_source_range = chain->reverse ? Interval(out.subject_range.end_, out.subject_range.begin_) : Interval(
                out.subject_range.begin_, out.subject_range.end_);
        out.frame = chain->reverse + 2;

        m.hsp.push_back(std::move(out));
    }

    return m;
}


// TODO: maybe include this into cigar_to_hsp; remove score computation from WfaCigar (makes more sense if short-mode is eventually removed); remove score from Cigar class
int compute_alignment_score(const Cigar &cigar, const Search::Config &cfg, const Sequence &target, const Sequence &query, const int pos_i, const int pos_j) {
    int score = 0;
    int pattern_pos = pos_i - cigar.max_query() - 1;
    int text_pos = pos_j - cigar.max_target() - 1;
    for (auto &c : cigar.cigar_data) {
        switch (c.second) {
            case 'M':
                // add score based on mismatch/match
                for (int i = 0; i < c.first; ++i) {
                    if (query[pattern_pos++] == target[text_pos++]) {
                        score += cfg.score_builder->reward();
                    } else {
                        score += cfg.score_builder->penalty();
                    }
                }
                break;
            case 'I':
                score -= cfg.score_builder->gap_open() + c.first * cfg.score_builder->gap_extend();
                pattern_pos += c.first;
                break;
            case 'D':
                score -= cfg.score_builder->gap_open() + c.first * cfg.score_builder->gap_extend();
                text_pos += c.first;
                break;
            default:
                break;
        }
    }


    return score;
}


Hsp build_align_HSP( const Search::Config &cfg, const BlockId id, const Chain *chain, const Sequence &query, const Sequence &target) {

    const auto kmer_size = shapes[0].length_;

    // extend to the left at start
    std::vector<Letter> query_left = query.subseq(0, chain->anchors.back().i).reverse();
    std::vector<Letter> target_left = target.subseq(std::max(0, (int) (chain->anchors.back().j - (query_left.size() * 2))),
                                                    chain->anchors.back().j).reverse();
    // left extension
    Cigar extension = config.dna_extension == DNAExtensionAlgo::WFA ?
                      static_cast<Cigar>(WfaCigar(target_left, query_left, cfg, true))
                                                                    :
                      KswCigar(target_left, query_left, cfg, KSW_FLAG_L);

    // TODO: reserve space for extension.cigar_data with + number of anchors?

    int anchor_distance_query= INT_MAX;
    int anchor_distance_target = INT_MAX;
    // iterate anchors from start to end (reverse)
    for (int i = static_cast<int>(chain->anchors.size()) - 1; i > 0; --i) {

        // if i no overlap with previous anchor, push match for anchor // TODO: maybe do it the other way around?
        if (anchor_distance_query > kmer_size && anchor_distance_target > kmer_size) {
            extension.cigar_data.emplace_back(kmer_size, 'M');
            //extension.add_score(kmer_size * cfg.score_builder->reward());
        }

        // compute distance of anchors
        anchor_distance_query = chain->anchors[i - 1].i - chain->anchors[i].i;
        anchor_distance_target = chain->anchors[i - 1].j - chain->anchors[i].j;

        // Case 1: Anchors overlap in Query and Target exactly the same
        if (anchor_distance_query == anchor_distance_target && anchor_distance_query <= kmer_size) {
            // push matches for last part of second anchor
            extension.cigar_data.emplace_back(anchor_distance_query, 'M');
            //extension.add_score(anchor_distance_query * cfg.score_builder->reward());
        }

            // Case 2: Anchors don't overlap in Query and Target
        else if (anchor_distance_query > kmer_size && anchor_distance_target > kmer_size) {
            // extend to the right between anchors
            Sequence query_right = query.subseq(chain->anchors[i].i + kmer_size, chain->anchors[i - 1].i);
            Sequence target_right = target.subseq(chain->anchors[i].j + kmer_size, chain->anchors[i - 1].j);

            extension = extension + (config.dna_extension == DNAExtensionAlgo::WFA ?
                                     static_cast<Cigar>(WfaCigar(target_right, query_right, cfg, false, true)) // takes a lot of memory, ask how to take less
                                                                                   :
                                     KswCigar(target_right, query_right, cfg, KSW_FLAG_B)); // TODO: maybe use adaptive KSW2_BAND? (based on anchor distance/difference in length of query and target)
        }


            // alternatively realign through these cases?

            // Case 3: Anchors overlap more in Query than in Target
        else if (anchor_distance_query <= kmer_size && anchor_distance_query < anchor_distance_target) {
            // push deletions
            int number_of_gaps = anchor_distance_target - anchor_distance_query;
            extension.cigar_data.emplace_back(number_of_gaps, 'D');
            //extension.add_score(- (number_of_gaps * cfg.score_builder->gap_extend()) - cfg.score_builder->gap_open());

            // push matches for last part of second anchor
            extension.cigar_data.emplace_back(anchor_distance_query, 'M');
            //extension.add_score(anchor_distance_query * cfg.score_builder->reward());
        }

            // Case 4: Anchors overlap more in Target than in Query
        else if (anchor_distance_target <= kmer_size && anchor_distance_target < anchor_distance_query) {
            // push insertions
            int number_of_gaps = anchor_distance_query - anchor_distance_target;
            extension.cigar_data.emplace_back(number_of_gaps, 'I');
            //extension.add_score(- (number_of_gaps * cfg.score_builder->gap_extend()) - cfg.score_builder->gap_open());

            // push matches for last part of second anchor
            extension.cigar_data.emplace_back(anchor_distance_target, 'M');
            //extension.add_score(anchor_distance_target * cfg.score_builder->reward());

        }
            // error message
        else {
            throw std::runtime_error("Error in chaining extension: No case matched");
        }
    } // end for loop

    // if no overlap with last anchor, push match for last anchor
    if (anchor_distance_query > kmer_size && anchor_distance_target > kmer_size) {
        extension.cigar_data.emplace_back(kmer_size, 'M');
        extension.add_score(kmer_size * cfg.score_builder->reward());
    }

    // extend to the right at end
    Sequence query_right = query.subseq(chain->anchors[0].i + kmer_size, query.length());
    Sequence target_right = target.subseq(chain->anchors[0].j + kmer_size, std::min(target.length(),
                                                                                         chain->anchors[0].j + kmer_size + (int) query_right.length() * 2));
    extension = extension + (config.dna_extension == DNAExtensionAlgo::WFA ?
                             static_cast<Cigar>(WfaCigar(target_right, query_right, cfg))
                                                                           :
                             KswCigar(target_right, query_right, cfg, KSW_FLAG_R));


    Hsp out = Hsp();
    //out.score = extension.score();
    out.score = compute_alignment_score(extension, cfg, target, query, chain->anchors.back().i, chain->anchors.back().j);
    out.bit_score = cfg.score_builder->blast_bit_Score(out.score);
    out.evalue = cfg.score_builder->blast_eValue(out.score, query.length());
    if (out.evalue >= config.max_evalue)
        return NULL;

    cigar_to_hsp(extension, cfg.target->seqs()[id], query, chain->anchors.back().i, chain->anchors.back().j, out, chain->reverse);

    return out;
}


/**
 * Extends a chain to an alignment
 * @param cfg
 * @param id
 * @param begin
 * @param end
 * @param query
 * @param map_score_threshold
 * @return match object
 */
Extension::Match target_extension_chaining(const Search::Config &cfg, const BlockId id, const Chain *begin, const Chain *end,
                                           const Sequence &query, const Sequence &query_reverse) {

    Extension::Match m = Extension::Match(id, cfg.target->seqs()[id], ::Stats::TargetMatrix(), 0, 0);

    const Sequence &target = cfg.target->seqs()[id];

    for (auto chain = begin; chain != end; ++chain) {

        m.hsp.push_back(chain->reverse ? build_align_HSP(cfg, id, chain, query_reverse, target)
                                       : build_align_HSP(cfg, id, chain, query, target));
    }

    return m;
}


/**
 * Computes the chains of a query
 * @param cfg
 * @param query
 * @param reverse
 * @param p
 * @return vector of chains
 */
std::vector<Chain> compute_chains(const Search::Config &cfg, const Sequence &query, bool reverse, const ChainingParameters &p) {
    std::vector<Chain> chains;


    auto hits = seed_lookup(query, cfg.target->seqs(), cfg.dna_ref_index.get(), cfg.minimizer_window);

    if (hits.empty()) return chains;

    // Sort SeedMatch vector by BlockId and target location j
    std::sort(hits.begin(), hits.end(), [](const SeedMatch &a, const SeedMatch &b) {
        return a.id() < b.id() || (a.id() == b.id() && a.j() < b.j());
    });


    auto hits_it = merge_keys(hits.begin(), hits.end(), [](const SeedMatch &hit1) { return hit1.id(); });

    // set minimum max distances
    // if (p.MAX_DIST_X < p.BAND_WIDTH) p.MAX_DIST_X = p.BAND_WIDTH;
    // if (p.MAX_DIST_Y < p.BAND_WIDTH) p.MAX_DIST_Y = p.BAND_WIDTH;

    // Iterate through the sorted hits vector and process subsets based on unique IDs
    while (hits_it.good()) {

        auto new_chains = chain_dp(cfg.minimizer_window, shapes[0].length_, p, hits_it.begin().operator->(),
                                   hits_it.end().operator->(), reverse);
        chains.insert(chains.end(), std::make_move_iterator(new_chains.begin()),
                      std::make_move_iterator(new_chains.end()));

        ++hits_it;
    }

    return chains;

}


/**
 * Chaining of a query and mapping/alignment of the chains
 * @param cfg
 * @param query
 * @return vector of matches
 */
std::vector<Extension::Match> query_extension_chaining(const Search::Config &cfg, const Sequence &query,  const Sequence &query_reverse) {
    std::vector<Extension::Match> matches;
    // cast cfg.chain_pen_gap to float
    ChainingParameters params(static_cast<float>(cfg.chain_pen_gap), static_cast<float>(cfg.chain_pen_skip));

    auto chains = compute_chains(cfg,query, false, params);
    auto chains_r = compute_chains(cfg,query_reverse, true, params);

    if (chains.empty() && chains_r.empty()) return matches;

    chains.insert( chains.end(), std::make_move_iterator(chains_r.begin()),
                   std::make_move_iterator(chains_r.end()));


    // sort chains by score
    std::sort(chains.begin(), chains.end(), std::greater<>());

    // compute chain quality
    compute_primary_chains(chains, shapes[0].length_);

    // filter chains by score for mapping // TODO: for alignment as well? keep all primary?
    if (config.chaining_out) {
        const int map_score_threshold = chains[0].chain_score * params.MAP_PERCENTAGE;
        auto lower_bound = std::lower_bound(chains.begin(), chains.end(), map_score_threshold, [](const Chain& chain, int threshold) {
            return chain.chain_score > threshold;
        });

        chains.erase(lower_bound, chains.end());
    }


    // sort chains by score and target id, will be necessary when primary for all targets
    std::sort(chains.begin(), chains.end(), [](const Chain &a, const Chain &b) {
        return (a.target_id == b.target_id) ? (a.chain_score > b.chain_score) : (a.target_id < b.target_id);
    });

    auto it_chains = merge_keys(chains.begin(), chains.end(),
                                [](const Chain &chain1) { return chain1.target_id; });

    while (it_chains.good()) {
        Extension::Match m = config.chaining_out ?
                             build_map_HSP(cfg, it_chains.key(), it_chains.begin().operator->(),
                                           it_chains.end().operator->()) :
                             target_extension_chaining(cfg, it_chains.key(), it_chains.begin().operator->(),
                                                       it_chains.end().operator->(), query, query_reverse);

        if (!m.hsp.empty()) {
            matches.push_back(std::move(m));
        }

        ++it_chains;
    }

    return matches;

}



}