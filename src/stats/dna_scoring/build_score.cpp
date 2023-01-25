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

#include <stdexcept>
#include "../lib/blast/blast_stat.h"
#include "../lib/blast/blast_encoding.h"
#include "../lib/blast/blast_setup.h"
#include "build_score.h"


namespace Stats {

Blastn_Score::Blastn_Score(const int reward, const int penalty, const int gapopen, const int gapextend,uint64_t db_letters,int64_t sequence_count) :
        m_ScoreBlk(BlastScoreBlkNew(BLASTNA_SEQ_CODE, 1)),
        reward_(reward),
        penalty_(penalty),
        gap_open_(gapopen),
        gap_extend_(gapextend),
        target_length_(db_letters),
        db_size_(sequence_count)

{
    int status;



    if(m_ScoreBlk == nullptr)
        throw std::runtime_error("Failed to initialize blast score block");


    m_ScoreBlk->kbp_gap_std[0] = Blast_KarlinBlkNew();
    Blast_ScoreBlkKbpIdealCalc(m_ScoreBlk);
    m_ScoreBlk->reward = reward_;
    m_ScoreBlk->penalty = penalty_;

    EBlastProgramType core_type =eBlastTypeBlastn;
    BlastScoringOptions *score_options;
    BlastScoringOptionsNew(core_type, &score_options);
    BLAST_FillScoringOptions(score_options, core_type, TRUE,penalty_,
                             reward_,
                             "",
                             gap_open_, gap_extend_);
    status = Blast_ScoreBlkMatrixInit(core_type, score_options,
                                      m_ScoreBlk, nullptr);
    score_options = BlastScoringOptionsFree(score_options);
    if (status)
        throw std::runtime_error("Failed to initialize scoring matrix");

    m_ScoreBlk->kbp_gap_std[0] = Blast_KarlinBlkNew();

    Blast_ScoreBlkKbpIdealCalc(m_ScoreBlk);
    status = Blast_KarlinBlkNuclGappedCalc(m_ScoreBlk->kbp_gap_std[0],
                                               gap_open_, gap_extend_,
                                               m_ScoreBlk->reward,
                                               m_ScoreBlk->penalty,
                                               m_ScoreBlk->kbp_ideal,
                                               &(m_ScoreBlk->round_down),
                                           nullptr);


    if (status || m_ScoreBlk->kbp_gap_std[0] == nullptr ||
    m_ScoreBlk->kbp_gap_std[0]->Lambda <= 0.0) {
        throw std::runtime_error("Failed to initialize Karlin Blocks");
    }

    kbp = m_ScoreBlk->kbp_gap_std[0] ;
}

double Blastn_Score::blast_bit_Score(int raw_score) const {

    return ((raw_score * kbp->Lambda) - kbp->logK) / NCBIMATH_LN2;
    }


double Blastn_Score::blast_eValue(int raw_score,int query_length) const {
    uint64_t searchspace = calculate_length_adjustment(query_length,query_length,this->target_length_) *calculate_length_adjustment(this->target_length_, query_length,this->target_length_,this->db_size_);
    return BLAST_KarlinStoE_simple(raw_score, kbp,(long)searchspace);
    }

uint64_t Blastn_Score::calculate_length_adjustment(uint64_t length, int query_length, uint64_t target_length, int64_t db_size) const {
    return length - (uint64_t)(expected_hsp_value(query_length,target_length) * (double)db_size);
}

uint64_t Blastn_Score::calculate_length_adjustment(uint64_t length, int query_length, uint64_t target_length) const {
    if(expected_hsp_value(query_length,target_length) < (1 / kbp->K))
        return (int64_t) (1/kbp->K);
    else
        return  length - (int64_t) expected_hsp_value(query_length,target_length);
}

double Blastn_Score::expected_hsp_value(int query_length, uint64_t target_length) const {
	return (log(kbp->K * query_length * (double)target_length) / kbp->H);
}

Blastn_Score::~Blastn_Score() {
	m_ScoreBlk = BlastScoreBlkFree(m_ScoreBlk);
}

}