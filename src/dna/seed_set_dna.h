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
#pragma once
namespace Dna{

struct SeedMatch{
    SeedMatch(Loc i, BlockId id, Loc j);
    int ungapped_score()const{return score_;}
    void l_ungapped(int l){this->l_ungapped_ = l;}
    void score(int r){this->score_ = r;}
    int l_ungapped()const{return l_ungapped_;}
    Loc i()const{return i_;}
    Loc j()const{return j_;}
    BlockId id()const{return target_id_;}
    int i_end() const {
        return i_ + score_ - 1;
    }
    int j_end() const {
        return j_ + score_ - 1;
    }
    void set_i(Loc new_i) {
        i_ = new_i;
    }
    void set_j(Loc new_j) {
        j_ = new_j;
    }

    bool operator>(const SeedMatch& hit)const
    {return this->id() < hit.id() || ((this->id() == hit.id()) && this->ungapped_score() > hit.ungapped_score());}


private:
    Loc i_;
    Loc j_;
    BlockId target_id_;
    int l_ungapped_;
    int score_;
};
std::vector<SeedMatch> seed_lookup(const Sequence &query, SequenceSet& target_seqs, const Index *filter, Loc it_param);
};
