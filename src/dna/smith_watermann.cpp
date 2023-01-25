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

#include "smith_watermann.h"
#include "../data/block/block.h"
#include "../data/sequence_file.h"
#include <iomanip>


using std::vector;


namespace SmithWaterman{

    int scoringFunction(int match, int mismatch, Letter first,Letter second){
        if(first == second)
            return match;
        else
            return mismatch;
    }

    int dynamic_programm(const Sequence &target, const Sequence &query, int match, int mismatch, int gapopen, int gapextend){
        vector<vector<int>> dp_matrix(query.length()+1,vector<int>(target.length()+1,0));
        vector<int> hgap(target.length()+1,0);

        int max = 0;
        for(int i = 1; i < (int)dp_matrix.size(); i++){
            int vgap = 0;
            for (int j = 1; j<(int)dp_matrix[0].size(); j++){
                int s = dp_matrix[i-1][j-1]+ scoringFunction(match,mismatch, query[i-1], target[j-1]);
                s = std::max({s,hgap[j],vgap,0},std::less<int>());
                dp_matrix[i][j] = s ;
                vgap -= gapextend;
                hgap[j] -= gapextend;
                int open = s - gapopen - gapextend;
                vgap = std::max(vgap, open);
                hgap[j] = std::max(hgap[j], open);

                if(dp_matrix[i][j] > max){
                    max = dp_matrix[i][j];

                }

            }
        }

        return max;
    }

    std::pair<vector<Extension::Match>, Extension::Stats> local_alignment(const Search::Config &cfg,BlockId query_id){
        vector<Extension::Match> matches;
        SequenceSet& target_seqs = cfg.target->seqs();
        Sequence query_sequence = cfg.query->seqs()[query_id];

        for(int i = 0; i < target_seqs.size(); i ++){
            Sequence target_sequence = target_seqs[i];
            Extension::Match m = Extension::Match(query_id,target_sequence,::Stats::TargetMatrix(),0,0);
            int score = dynamic_programm(target_sequence, query_sequence,cfg.score_builder->reward(), cfg.score_builder->penalty(), cfg.score_builder->gap_open(), cfg.score_builder->gap_extend()) ;
            m.hsp.emplace_back(Hsp(false,score));
            m.hsp.back().bit_score = cfg.score_builder->blast_bit_Score(score);
            m.hsp.back().evalue = cfg.score_builder->blast_eValue(score,query_sequence.length());
            matches.emplace_back(m);

        }
        return { matches, Extension::Stats()};
    }

}