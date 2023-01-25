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

#include "blast/blast_stat.h"

namespace Stats {
struct Blastn_Score{
    Blastn_Score(const int reward,const int penalty,const  int gapopen,const int gapextend,uint64_t db_letters,int64_t sequence_count);
    double blast_bit_Score(int raw_score) const;
    double blast_eValue(int raw_score,int query_length) const;
    ~Blastn_Score();
    int reward()const{return reward_;}
    int penalty()const{return penalty_;}
    int gap_open()const{return gap_open_;}
    int gap_extend()const{return gap_extend_;}


private:
    Blast_KarlinBlk* kbp;
    BlastScoreBlk *m_ScoreBlk;
    int reward_,penalty_,gap_open_,gap_extend_;
    uint64_t target_length_;
    int64_t db_size_;
    uint64_t calculate_length_adjustment(uint64_t length,int query_length, uint64_t target_length, int64_t db_size) const;
    uint64_t calculate_length_adjustment(uint64_t length,int query_length, uint64_t target_length) const;
    double expected_hsp_value(int query_length, uint64_t target_length) const;
};
}
