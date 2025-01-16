/****
DIAMOND protein aligner
Copyright (C) 2024 Vincent Spath

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

#include "extension_seed_matches.h"
#include "../data/block/block.h"
#include "../util/util.h"

namespace Dna {

std::vector<SeedMatch> extend_seeds_ungapped(const SeedMatch *seedMatchBegin, const SeedMatch *seedMatchEnd, const Sequence &query, const Search::Config &cfg) {
    const Sequence &target = cfg.target->seqs()[seedMatchBegin->id()];
    const int kmer_size = shapes[0].length_;
    std::vector<SeedMatch> new_hits;
    new_hits.reserve(std::distance(seedMatchBegin, seedMatchEnd));

    int prev_diagonal = std::numeric_limits<int>::min();

    const int query_length = query.length();
    const int target_length = target.length();
    for (auto hit = seedMatchBegin; hit != seedMatchEnd; ++hit) {

        const int diagonal = hit->i() - hit->j();

        if (diagonal == prev_diagonal && hit->i() <= new_hits.back().i()) {
            continue;
        }

        int left_i = hit->i() - 1;
        int left_j = hit->j() - 1;
        while (query[left_i] == target[left_j] && left_i > -1 && left_j > -1) {
            left_i--;
            left_j--;
        }

        int right_i = hit->i() + kmer_size;
        int right_j = hit->j() + kmer_size;
        while (query[right_i] == target[right_j] && right_i < query_length && right_j < target_length) {
            right_i++;
            right_j++;
        }

        new_hits.emplace_back(right_i, hit->id(), right_j);
        new_hits.back().score(right_i - left_i - 1);

        prev_diagonal = diagonal;
    }

    return new_hits;
}


bool compare_by_target_and_diagonal(const SeedMatch& a, const SeedMatch& b) {
    if (a.id() != b.id()) {
        return a.id() < b.id();
    }
    if ((a.i() - a.j()) != (b.i() - b.j())) {
        return (a.i() - a.j()) < (b.i() - b.j());
    }
    return a.i() < b.i();
}


std::vector<SeedMatch>
merge_and_extend_seeds(std::vector<SeedMatch> &seed_hits, const Sequence &query, const Search::Config &cfg) {

    std::sort(seed_hits.begin(), seed_hits.end(), compare_by_target_and_diagonal);

    auto hits_it_target = merge_keys(seed_hits.begin(), seed_hits.end(), [](const SeedMatch &hit) { return hit.id(); });

    std::vector<SeedMatch> new_seed_hits;
    new_seed_hits.reserve(seed_hits.size());

    while (hits_it_target.good()) {
        auto new_hits = extend_seeds_ungapped(hits_it_target.begin().operator->(), hits_it_target.end().operator->(), query, cfg);

        new_seed_hits.insert(new_seed_hits.end(), std::make_move_iterator(new_hits.begin()), std::make_move_iterator(new_hits.end())); //move sinnvoll?
        ++hits_it_target;
    }

    return new_seed_hits;
}

}
