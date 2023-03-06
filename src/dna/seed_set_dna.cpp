/****
DIAMOND protein aligner
Copyright (C) 2022 Dimitrios Koutsogiannis

Code developed by Dimitrios Koutsogiannis <dimitrios.koutsogiannis@tue.mpg.de>

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
#include <utility>
#include "../basic/shape_config.h"
#include "../basic/seed_iterator.h"
#include "../data/block/block.h"
#include "../data/sequence_file.h"

using std::pair;

namespace Dna {

std::vector<SeedMatch> seed_lookup(const Sequence &query, SequenceSet& target_seqs, const Index *filter, Loc window_size){
    std::vector<SeedMatch> out;
    std::vector<Letter> buf = query.copy();

    const Shape &sh = shapes[0];
    MinimizerIterator it(buf, sh, window_size);

    while (it.good()) {
        uint64_t key = *it;
        pair<SeedArray::Entry*,SeedArray::Entry*> ref_position = filter->contains(key);
        if (ref_position.first != nullptr){
            for(auto pos = ref_position.first; pos != ref_position.second; pos++){
                pair<BlockId ,Loc> id_loc = target_seqs.local_position(pos->value);
                out.emplace_back(it.pos(), id_loc.first, id_loc.second); //ref_position));
            }
        }
        ++it;
    }

    return out;
}


SeedMatch::SeedMatch(Loc i, BlockId id, Loc j):
i_(i),
target_id_(id),
j_(j)
{}
}
