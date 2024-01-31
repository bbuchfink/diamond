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

#include "../lib/ksw2/ksw2.h"
#include "seed_set_dna.h"
#include "../util/util.h"
#include "bindings/cpp/WFAligner.hpp"
#include "extension.h"
#include "timer.h"
#include "chain.h"
#include "extension_chain.h"
#include "../basic/shape_config.h"
#include "../align/target.h"

using std::make_move_iterator;

const EMap<DNAExtensionAlgo> EnumTraits<DNAExtensionAlgo>::to_string = { {DNAExtensionAlgo::KSW,"ksw"}, { DNAExtensionAlgo::WFA,"wfa"} };
const SEMap<DNAExtensionAlgo> EnumTraits<DNAExtensionAlgo>::from_string = { {"ksw",DNAExtensionAlgo::KSW}, {"wfa",DNAExtensionAlgo::WFA} };

namespace Dna {


    struct ExtendedSeed {
        ExtendedSeed(int i_min, int i_max, int j_min, int j_max) :
                i_min_(i_min),
                i_max_(i_max),
                j_min_(j_min),
                j_max_(j_max) {}

        int i_min_extended() const { return i_min_; };

        int i_max_extended() const { return i_max_; };

        int j_min_extended() const { return j_min_; };

        int j_max_extended() const { return j_max_; };

    private:
        int i_min_, i_max_, j_min_, j_max_;
    };



    bool intersection(SeedMatch &hit, const std::vector<ExtendedSeed> &extended, int window) {
        return std::any_of(extended.begin(), extended.end(), [hit, window](ExtendedSeed s) {
            return hit.i() + 1 >= s.i_min_extended() && hit.i() + window <= s.i_max_extended() &&
                   hit.j() + 1 >= s.j_min_extended() && hit.j() + window <= s.j_max_extended();
        });
    }

    void calculate_ungapped_scores(SeedMatch &hit, const Sequence &target, const Sequence &query) {
        int score = 0;
        int i = 0;

        while (query[hit.i() - i] == target[hit.j() - i]) {
            score++;
            i++;
        }
        hit.l_ungapped(score);

        i = 1;
        while (query[hit.i() + i] == target[hit.j() + i]) {
            score++;
            i++;
        }
        hit.score(score);

    }


    KswCigar::KswCigar(const Sequence &tseq, const Sequence &qseq, const Search::Config &cfg, int flag) {


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
        ksw_extz(nullptr, qseq.length(), reinterpret_cast<const uint8_t *>(qseq.data()), tseq.length(),
                 reinterpret_cast<const uint8_t *>(tseq.data()), NUCLEOTIDE_COUNT, mat,
                 (int8_t) cfg.score_builder->gap_open(), (int8_t) cfg.score_builder->gap_extend(), KSW2_BAND, //config.padding
                 config.zdrop, flag, &ez);
#else

        ksw_extz2_sse(nullptr, qseq.length(), reinterpret_cast<const uint8_t *>(qseq.data()), tseq.length(),
                  reinterpret_cast<const uint8_t *>(tseq.data()), NUCLEOTIDE_COUNT, mat,
                  (int8_t) cfg.score_builder->gap_open(), (int8_t) cfg.score_builder->gap_extend(), KSW2_BAND,
                  config.zdrop, KSW2_END_BONUS, flag, &ez);

#endif
        this->score_ = ez.max;
        this->max_query_ = ez.max_q;
        this->max_target_ = ez.max_t;
        cigar_data.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar_data.emplace_back(ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);


        free(ez.cigar);

    }



    WfaCigar::WfaCigar(const Sequence &tseq, const Sequence &qseq, const Search::Config &cfg, bool left, bool global) {

        std::string tseq2 = tseq.to_string();
        std::string qseq2 = qseq.to_string();

        int tl = tseq2.length(), ql = qseq2.length();

        thread_local std::unique_ptr<wfa::WFAlignerGapAffine> aligner;
        if (!aligner)
            aligner.reset(new wfa::WFAlignerGapAffine(0, -config.mismatch_penalty, config.gap_open, config.gap_extend,
                                                      wfa::WFAligner::Alignment, wfa::WFAligner::MemoryLow)); // TODO: try memory modes / use dependent on length
        aligner->setHeuristicNone();
        aligner->setHeuristicZDrop(config.zdrop, WFA_CUTOFF_STEPS);
        aligner->setHeuristicWFadaptive(10,50,1);

        if (global) {
            aligner->alignEnd2End(tseq2.c_str(), tl, qseq2.c_str(), ql);
            //aligner->alignEnd2End(tseq2, qseq2);  // TODO: global alignment does not work properly sometimes (update wfa-lib when possible)
        } else {
            aligner->alignExtension(tseq2.c_str(), tl, qseq2.c_str(), ql);
        }

        auto cigar = aligner->getCIGAR(true);


        max_query_ = -1;
        max_target_ = -1;
        int steps = 0;
        for (char c: cigar) {
            if (isdigit(c)) {
                steps = steps * 10 + (c - '0');
                continue;
            }
            cigar_data.emplace_back(steps, c);
            switch (c) {
                case '=':
                    this->score_ += cfg.score_builder->reward() * steps;
                    max_query_ += steps;
                    max_target_ += steps;
                    break;
                case 'X':
                    this->score_ += cfg.score_builder->penalty() * steps;
                    max_query_ += steps;
                    max_target_ += steps;
                    break;
                case 'I':
                    this->score_ -= cfg.score_builder->gap_open() + steps * cfg.score_builder->gap_extend();
                    max_query_ += steps;
                    break;
                case 'D':
                    this->score_ -= cfg.score_builder->gap_open() + steps * cfg.score_builder->gap_extend();
                    max_target_ += steps;
                    break;
                default:
                    throw std::runtime_error(std::string("WFA Cigar: Invalid Cigar Symbol ") + c);

            }
            steps = 0;
        }
        if (left)
            std::reverse(cigar_data.begin(), cigar_data.end());
    }



    void cigar_to_hsp(const Sequence &target, const Sequence &query, const SeedMatch &hit, Hsp &out, bool reverse) {
        int pattern_pos = hit.i() - hit.l_ungapped() +1;
        int text_pos = hit.j() - hit.l_ungapped() +1;
        out.query_range.begin_ = pattern_pos;
        out.subject_range.begin_ = text_pos;

        for (int i = 0; i < hit.ungapped_score(); ++i) {
            out.push_match(target[text_pos], query[pattern_pos], true);
            pattern_pos++;
            text_pos++;
        }

        out.query_range.end_ = pattern_pos;
        out.subject_range.end_ = text_pos;
        out.transcript.push_terminator();
        out.target_seq = target;
        out.query_source_range = out.query_range;
        out.subject_source_range = reverse ? Interval(out.subject_range.end_, out.subject_range.begin_) : Interval(
                out.subject_range.begin_, out.subject_range.end_);
        out.frame = reverse + 2;

    }


    void cigar_to_hsp(const Cigar &cigar, const Sequence &target, const Sequence &query, const int pos_i, const int pos_j, Hsp &out,
                      bool reverse) {

        int pattern_pos = pos_i - cigar.max_query() - 1;
        int text_pos = pos_j - cigar.max_target() - 1;
        out.query_range.begin_ = pattern_pos;
        out.subject_range.begin_ = text_pos;

        for (auto operation: cigar.cigar_data) {
            switch (operation.second) {
                case 'M':
                case '=':
                case 'X':
                    for (int j = 0; j < operation.first; ++j) {
                        out.push_match(target[text_pos], query[pattern_pos], true);
                        pattern_pos++;
                        text_pos++;
                    }
                    break;
                case 'D':
                    out.push_gap(op_deletion, operation.first, target.data() + operation.first + text_pos);
                    text_pos += operation.first;
                    break;
                case 'I':
                    out.push_gap(op_insertion, operation.first, query.data() + operation.first + pattern_pos);
                    pattern_pos += operation.first;
                    break;
                default:
                    break;
            }
        }

        out.query_range.end_ = pattern_pos;
        out.subject_range.end_ = text_pos;
        out.transcript.push_terminator();
        out.target_seq = target;
        out.query_source_range = out.query_range;
        out.subject_source_range = reverse ? Interval(out.subject_range.end_, out.subject_range.begin_) : Interval(
                out.subject_range.begin_, out.subject_range.end_);
        out.frame = reverse + 2;
    }



    Extension::Match
    target_extension(const Search::Config &cfg, const BlockId id, const Sequence &query, SeedMatch *begin,
                     SeedMatch *end, bool reverse) {
        std::vector<ExtendedSeed> extended{};

        Extension::Match m = Extension::Match(id, cfg.target->seqs()[id], ::Stats::TargetMatrix(), 0, 0);
        for (auto hit = begin; hit != end; ++hit) {
            if ((intersection(*hit, extended, cfg.minimizer_window)))
                continue;
            const Sequence &target = cfg.target->seqs()[id];
            if (hit->ungapped_score() == query.length()) {
                Hsp out = Hsp();

                out.score = hit->ungapped_score() * cfg.score_builder->reward();
                out.bit_score = cfg.score_builder->blast_bit_Score(out.score);
                out.evalue = cfg.score_builder->blast_eValue(out.score, query.length());
                if (out.evalue >= config.max_evalue)
                    continue;

                cigar_to_hsp(cfg.target->seqs()[id], query, *hit, out, reverse);
                m.hsp.push_back(out);
                extended.emplace_back(hit->i() - hit->l_ungapped(),
                                      out.query_range.end_,
                                      hit->j() - hit->l_ungapped(),
                                      out.subject_range.end_);

            } else {

                Sequence query_right = query.subseq(hit->i(), query.length());
                Sequence target_right = target.subseq(hit->j(), std::min(target.length(),
                                                                         hit->j() + (int) query_right.length() * 2));

                std::vector<Letter> query_left = query.subseq(0, hit->i()).reverse();
                std::vector<Letter> target_left = target.subseq(std::max(0, (int) (hit->j() - (query_left.size() * 2))),
                                                                hit->j()).reverse();


                Cigar extension = config.dna_extension == DNAExtensionAlgo::WFA ?
                                  static_cast<Cigar>(WfaCigar(target_left, query_left, cfg, 1) +
                                                     WfaCigar(target_right, query_right, cfg))
                                                                                :
                                  KswCigar(target_left, query_left, cfg, KSW_FLAG_L) +
                                  KswCigar(target_right, query_right, cfg, KSW_FLAG_R);

                Hsp out = Hsp();
                out.score = extension.score();
                out.bit_score = cfg.score_builder->blast_bit_Score(out.score);
                out.evalue = cfg.score_builder->blast_eValue(out.score, query.length());
                if (out.evalue >= config.max_evalue)
                    continue;

                cigar_to_hsp(extension, cfg.target->seqs()[id], query, hit->i(),hit->j(), out, reverse);

                m.hsp.push_back(out);
                extended.emplace_back(hit->i() - extension.max_query(),
                                      out.query_range.end_,
                                      hit->j() - extension.max_target(),
                                      out.subject_range.end_);

            }
        }
        m.inner_culling();
        return m;
    }


    std::vector<Extension::Match> query_extension(const Search::Config &cfg, const Sequence &query, bool reverse) {
        std::vector<Extension::Match> matches;


        auto hits = seed_lookup(query, cfg.target->seqs(), cfg.dna_ref_index.get(), cfg.minimizer_window);

        // calculate ungapped scores for all hits
        std::for_each(hits.begin(), hits.end(), [&](SeedMatch &hit) {
            calculate_ungapped_scores(hit, cfg.target->seqs()[hit.id()], query);
        });

        std::sort(hits.begin(), hits.end(), std::greater<>());

        auto it = merge_keys(hits.begin(), hits.end(), [](const SeedMatch &hit1) { return hit1.id(); });

        while (it.good()) {
            Extension::Match m = target_extension(cfg, it.key(), query, it.begin().operator->(), it.end().operator->(),
                                                  reverse);
            if (!m.hsp.empty()) {
                matches.push_back(m);
            }

            ++it;
        }

        return matches;
    }





        std::pair<std::vector<Extension::Match>, Extension::Stats> extend(const Search::Config &cfg,const Sequence &query) {

        // perform chaining (only mapping atm)
        if (config.chaining_out || config.align_long_reads) {
           return {query_extension_chaining(cfg,query, Translator::reverse(query)), Extension::Stats()};
        }
        else {

        std::vector<Extension::Match> matches = query_extension(cfg,query, false);
        std::vector<Extension::Match> reverse = query_extension(cfg, Translator::reverse(query), true);

        matches.insert( matches.end(), make_move_iterator(reverse.begin()), make_move_iterator(reverse.end()));

        culling(matches, cfg);

        return {matches, Extension::Stats()};
        }

    }
}






