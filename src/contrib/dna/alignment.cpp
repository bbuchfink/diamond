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

#include "alignment.h"
#include "../lib/ksw2/ksw2.h"
#include "seed_set_dna.h"
#include "../util/util.h"
#include "bindings/cpp/WFAligner.hpp"
#include "timer.h"
#include "chain.h"
#include "extension_chain.h"
#include "../align/target.h"



namespace Dna {


AlignmentStatus
compute_ksw_cigar(const Sequence &targetSequence, const Sequence &querySequence, const Search::Config &cfg, int flag,
                  Cigar &extension, const int zdrop, const int band) {

    int a = cfg.score_builder->reward(), b = cfg.score_builder->penalty() < 0 ? cfg.score_builder->penalty()
                                                                              : -cfg.score_builder->penalty(); // a>0 and b<0
    int8_t mat[NUCLEOTIDE_COUNT * NUCLEOTIDE_COUNT] = {static_cast<int8_t>(a), static_cast<int8_t>(b),
                                                       static_cast<int8_t>(b),
                                                       static_cast<int8_t>(b), 0, static_cast<int8_t>(b),
                                                       static_cast<int8_t>(a),
                                                       static_cast<int8_t>(b), static_cast<int8_t>(b), 0,
                                                       static_cast<int8_t>(b),
                                                       static_cast<int8_t>(b), static_cast<int8_t>(a),
                                                       static_cast<int8_t>(b), 0,
                                                       static_cast<int8_t>(b), static_cast<int8_t>(b),
                                                       static_cast<int8_t>(b),
                                                       static_cast<int8_t>(a), 0, 0, 0, 0, 0, 0};
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
#ifdef __APPLE__
    ksw_extz(nullptr, querySequence.length(), reinterpret_cast<const uint8_t *>(querySequence.data()), targetSequence.length(),
             reinterpret_cast<const uint8_t *>(targetSequence.data()), NUCLEOTIDE_COUNT, mat,
             (int8_t) cfg.score_builder->gap_open(), (int8_t) cfg.score_builder->gap_extend(), band, //config.padding
             zdrop, flag, &ez);
#else
    ksw_extz2_sse(nullptr, querySequence.length(), reinterpret_cast<const uint8_t *>(querySequence.data()), targetSequence.length(),
              reinterpret_cast<const uint8_t *>(targetSequence.data()), NUCLEOTIDE_COUNT, mat,
              (int8_t) cfg.score_builder->gap_open(), (int8_t) cfg.score_builder->gap_extend(), band,
              zdrop, KSW2_END_BONUS, flag, &ez);
#endif

    // ez.max is extension (peak) score, ez.score is global alignment score
    const auto alignment_status = static_cast<AlignmentStatus>(ez.zdropped);

    if (flag == KSW_FLAG_L || flag == KSW_FLAG_R || alignment_status == DROPPED) {
        extension.score += ez.max;
    } else {
        extension.score += ez.score;
    }


    if (flag == KSW_FLAG_L)
        extension.setMaxValues(ez.max_q, ez.max_t);
    else if (extension.score < 1) {
        return NEGATIVE_SCORE;
    }

    for (int i = 0; i < ez.n_cigar; ++i)
        extension.extendCigar(ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);


    free(ez.cigar);

    return alignment_status;
}

std::pair<AlignmentStatus, std::string>
compute_wfa_extension(const std::string &querySequence, const std::string &targetSequence, const int band) {
    thread_local std::unique_ptr<wfa::WFAlignerGapAffine> aligner_extension;
    if (!aligner_extension) {
        aligner_extension.reset(new wfa::WFAlignerGapAffine(0, -config.mismatch_penalty, config.gap_open, config.gap_extend,
                                                  wfa::WFAligner::Alignment));

        aligner_extension->setHeuristicNone();
        aligner_extension->setHeuristicWFadaptive(10,50,1);
        //aligner_extension->setHeuristicZDrop(WFA_ZDROP_EXTENSION, WFA_CUTOFF_STEPS);
        aligner_extension->setHeuristicXDrop(100, 1);

        //aligner_extension->setHeuristicBandedAdaptive(-band,band,1);
    }

    aligner_extension->alignExtension(targetSequence.c_str(), targetSequence.length(), querySequence.c_str(), querySequence.length());

    return std::make_pair(static_cast<AlignmentStatus>(aligner_extension->getAlignmentStatus()),
                          aligner_extension->getCIGAR(true));
}




std::pair<AlignmentStatus, std::string>
compute_wfa_global(const std::string &querySequence, const std::string &targetSequence, int band) {

    band = std::min(band, WFA_BAND_EXTENSION);

    thread_local std::unique_ptr<wfa::WFAlignerGapAffine> aligner_global;
    if (!aligner_global) {
        aligner_global.reset(new wfa::WFAlignerGapAffine(0, -config.mismatch_penalty, config.gap_open, config.gap_extend,
                                                         wfa::WFAligner::Alignment));

        aligner_global->setHeuristicNone();
        aligner_global->setHeuristicWFadaptive(10,50,1);
        //aligner_global->setHeuristicZDrop(WFA_ZDROP_GLOBAL, WFA_CUTOFF_STEPS);
        aligner_global->setHeuristicXDrop(100, 1);
        //aligner_global->setHeuristicBandedAdaptive(-band,band,1);
    }



    aligner_global->alignEnd2End(targetSequence.c_str(), targetSequence.length(), querySequence.c_str(), querySequence.length());

    return std::make_pair(static_cast<AlignmentStatus>(aligner_global->getAlignmentStatus()),
                          aligner_global->getCIGAR(true));
}

AlignmentStatus
compute_wfa_cigar(const Search::Config &cfg, const std::string &querySequence, Cigar &extension, bool left, bool global,
                  const std::string &targetSequence, const int band) {

    auto [alignment_status, cigar] = global ?
            compute_wfa_global(querySequence, targetSequence, band)
            :
            compute_wfa_extension(querySequence, targetSequence, band);



    std::vector<std::pair<int, char>> cigar_data;
    // TODO: gibt es eine effizientere (aber elegante) Methode, als jedes Mal max_query und max_target zu berechnen?
    int max_query = -1;
    int max_target = -1;
    int steps = 0;
    for (char c: cigar) {
        if (isdigit(c)) {
            steps = steps * 10 + (c - '0');
            continue;
        }
        cigar_data.emplace_back(steps, c);
        switch (c) {
            case '=':
                extension.score += steps * cfg.score_builder->reward();
                max_query += steps;
                max_target += steps;
                break;
            case 'X':
                extension.score += steps * cfg.score_builder->penalty();
                max_query += steps;
                max_target += steps;
                break;
            case 'I':
                extension.score -= cfg.score_builder->gap_open() + (steps * cfg.score_builder->gap_extend());
                max_query += steps;
                break;
            case 'D':
                extension.score -= cfg.score_builder->gap_open() + (steps * cfg.score_builder->gap_extend());
                max_target += steps;
                break;
            default:
                throw std::runtime_error(std::string("WFA Cigar_short: Invalid Cigar_short Symbol ") + c);

        }

        steps = 0;
    }

    if (left) {
        std::reverse(cigar_data.begin(), cigar_data.end());
        extension.setMaxValues(max_query, max_target);
    }
    else if (extension.score < 1) {
        return NEGATIVE_SCORE;
    }
    extension.extendCigar(cigar_data);

    // return z-dropped
    return alignment_status;
}


Hsp build_hsp_from_cigar(const Cigar &cigar, const Sequence &target, const Sequence &query, const int firstAnchor_i,
                         const int firstAnchor_j, bool isReverse, const Search::Config &cfg) {
// timer
    ExtensionTimer timer_build;
    auto start_build = std::chrono::high_resolution_clock::now();

    Hsp alignHsp = Hsp();

    int query_pos = firstAnchor_i - cigar.query_extension_distance() - 1;
    int target_pos = firstAnchor_j - cigar.target_extension_distance() - 1;
    alignHsp.query_range.begin_ = query_pos;
    alignHsp.subject_range.begin_ = target_pos;

    for (auto& operation: cigar.getCigarDataConst()) {
        switch (operation.second) {
            case 'M':
            case '=':
            case 'X':
                for (int j = 0; j < operation.first; ++j) {
                    alignHsp.push_match(target[target_pos++], query[query_pos++], true);
                }
                break;
            case 'D':
                alignHsp.push_gap(op_deletion, operation.first, target.data() + operation.first + target_pos);
                target_pos += operation.first;
                break;
            case 'I':
                alignHsp.push_gap(op_insertion, operation.first, query.data() + operation.first + query_pos);
                query_pos += operation.first;
                break;
            default:
                break;
        }
    }

    alignHsp.score = cigar.score;
    alignHsp.bit_score = cfg.score_builder->blast_bit_Score(alignHsp.score);
    alignHsp.evalue = cfg.score_builder->blast_eValue(alignHsp.score, query.length());
    if (alignHsp.evalue >= config.max_evalue)
                return alignHsp;
    alignHsp.query_range.end_ = query_pos;
    alignHsp.subject_range.end_ = target_pos;
    alignHsp.transcript.push_terminator();
    alignHsp.target_seq = target;
    alignHsp.query_source_range = alignHsp.query_range;
    alignHsp.subject_source_range = isReverse ? Interval(alignHsp.subject_range.end_, alignHsp.subject_range.begin_) : Interval(
            alignHsp.subject_range.begin_, alignHsp.subject_range.end_);
    alignHsp.frame = isReverse;//+ 2;
//time
    auto end_build = std::chrono::high_resolution_clock::now();
    timer_build.update(2,end_build-start_build);
    {
        std::lock_guard<std::mutex> lock(cfg.timer->mtx);
        *(cfg.timer) += timer_build;
    }

    return alignHsp;
}

}