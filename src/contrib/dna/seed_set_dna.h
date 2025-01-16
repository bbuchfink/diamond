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

#include "../data/flags.h"
#include "dna_index.h"
#include "../basic/shape_config.h"
#pragma once
namespace Dna{

struct SeedMatch{
    SeedMatch(Loc i, BlockId id, Loc j);
    int ungapped_score()const{return score_;}
    void score(int scr){this->score_ = scr;}
    Loc i()const{return i_;}
    Loc j()const{return j_;}
    BlockId id()const{return target_id_;}
    int i_start() const {
        return i_ - score_;
    }
    int j_start() const {
        return j_ - score_;
    }

    bool operator>(const SeedMatch& hit)const
    {return this->id() < hit.id() || ((this->id() == hit.id()) && this->ungapped_score() > hit.ungapped_score());}


private:
    Loc i_;
    Loc j_;
    BlockId target_id_;
    int score_ = shapes[0].length_;
};
std::vector<SeedMatch> seed_lookup(const Sequence &query, SequenceSet& target_seqs, const Index *filter, Loc it_param);
};
