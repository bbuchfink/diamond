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

#include <string.h>
#include "../lib/ksw2/ksw2.h"
#include "ksw2_extension.h"
#include "seed_set_dna.h"
#include "../util/util.h"
#include "bindings/cpp/WFAligner.hpp"


namespace Dna {
struct ExtendedSeed{
    ExtendedSeed(int i_min, int i_max, int j_min, int j_max):
            i_min_(i_min),
            i_max_(i_max),
            j_min_(j_min),
            j_max_(j_max)
    {}

    int i_min_extended()const{return i_min_;};
    int i_max_extended()const{return i_max_;};
    int j_min_extended()const{return j_min_;};
    int j_max_extended()const{return j_max_;};

private:
    int i_min_, i_max_, j_min_, j_max_;
};

/**
ExtensionTimer::ExtensionTimer():
total_time(),
seed_lookup_time(),
ungapped(),
trimming(),
ksw_extension(),
diamond_cigar()
{
}

void ExtensionTimer::update(int operation, std::chrono::duration<long, std::ratio<1, 1000000000>> duration) {
    switch(operation){
        case 0:
            total_time.operator+=(duration);
            break;
        case 1:
            seed_lookup_time.operator+=(duration);
            break;
        case 2:
            ungapped.operator+=(duration);
            break;
        case 3:
            trimming.operator+=(duration);
            break;
        case 4:
            ksw_extension.operator+=(duration);
            break;
        case 5:
            diamond_cigar.operator+=(duration);
            break;
    }
}
 **/

struct KswCigar{
    KswCigar(ksw_extz_t ez):
        score_(ez.max),
        max_query_(ez.max_q),
        max_target_(ez.max_t)
        {
        cigar_data.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar_data.emplace_back(ez.cigar[i]>>4, "MID"[ez.cigar[i] & 0xf]);
        }

    KswCigar& operator+(const KswCigar& other){
       this->cigar_data.insert(this->cigar_data.end(), other.cigar_data.begin(), other.cigar_data.end());
       this->score_ += other.score_;
       return *this;
    }
    int32_t score()const{return score_;}
    int max_query()const{return max_query_;}
    int max_target()const{return max_target_;}
    std::vector<std::pair<int,char>> cigar_data;

private:
    int32_t score_;
    const int max_query_;
    const int max_target_;
};

bool intersection(SeedMatch &hit, const std::vector<ExtendedSeed> &extended){
    return std::any_of(extended.begin(), extended.end(), [hit](ExtendedSeed s)
    {return hit.i() >= s.i_min_extended() && hit.i() + 15 <= s.i_max_extended() && hit.j() >= s.j_min_extended() && hit.j()+15 <= s.j_max_extended() ;});
}

void align_wfa(const Sequence& tseq, const Sequence& qseq, int sc_mch, int sc_mis, int gapo, int gape) {
    // Parameters
    std::string tseq2 = tseq.to_string();
    std::string qseq2 = qseq.to_string();


    int tl = tseq2.length(), ql = qseq2.length();

    // Create aligner
    wfa::WFAlignerGapAffine aligner(-sc_mch,sc_mis,gapo,gape,wfa::WFAligner::Alignment);
    aligner.setHeuristicNone();
    aligner.setHeuristicZDrop(40,1); //20

    // Align sequences
    aligner.alignExtension(tseq2.c_str(),tl,qseq2.c_str(),ql);

    std::cout<<"WFA: " <<aligner.getCIGARString(false) << "\n";

}

KswCigar ksw2_align(const Sequence& tseq, const Sequence& qseq, int sc_mch, int sc_mis, int gapo, int gape,int flag) {
    int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[NUCLEOTIDE_COUNT * NUCLEOTIDE_COUNT] = {static_cast<int8_t>(a), static_cast<int8_t>(b), static_cast<int8_t>(b),
                      static_cast<int8_t>(b), 0, static_cast<int8_t>(b), static_cast<int8_t>(a),
                      static_cast<int8_t>(b), static_cast<int8_t>(b), 0, static_cast<int8_t>(b),
                      static_cast<int8_t>(b), static_cast<int8_t>(a), static_cast<int8_t>(b), 0,
                      static_cast<int8_t>(b), static_cast<int8_t>(b), static_cast<int8_t>(b),
                      static_cast<int8_t>(a), 0, 0, 0, 0, 0, 0};
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    ksw_extz2_sse(nullptr, qseq.length(), reinterpret_cast<const uint8_t *>(qseq.data()), tseq.length(),reinterpret_cast<const uint8_t *>(tseq.data()), NUCLEOTIDE_COUNT, mat, (int8_t)gapo, (int8_t)gape, config.padding == 0 ? -1 : config.padding , config.zdrop , KSW2_END_BONUS, flag, &ez);
    KswCigar out(ez);
    free(ez.cigar);


    return out;
}
void calculate_ungapped_scores(SeedMatch &hit, const Sequence &target, const Sequence &query){
    int score = 0;
    int i = 0;

    while(query[hit.i() - i] == target[hit.j()-i]){
        score++;
        i++;
    }

    i = 1 ;
    while(query[hit.i() + i] == target[hit.j()+i]){
        score++;
        i++;
    }
    hit.ungapped_score(score);
}
void cigar_to_hsp(const KswCigar& cigar, const Sequence &target, const Sequence &query, const SeedMatch& hit, Hsp &out) {

    int pattern_pos = hit.i() - cigar.max_query() -1;
    int text_pos = hit.j() - cigar.max_target() -1 ;
    out.query_range.begin_ = hit.i() - cigar.max_query() -1;
    out.subject_range.begin_ = hit.j() - cigar.max_target() -1;

    for(auto operation: cigar.cigar_data){
        switch(operation.second){
            case 'M':
                for(int j = 0; j < operation.first; ++j){
                    out.push_match(target[text_pos], query[pattern_pos],true);
                    pattern_pos++; text_pos++;
                }
                break;
            case 'D':
                out.push_gap(op_deletion,operation.first,target.data()+operation.first+text_pos);
                text_pos += operation.first;
                break;
            case 'I':
                out.transcript.push_back(op_insertion, (unsigned)operation.first);
                pattern_pos += operation.first;
                break;
            default:
                break;
        }

    }


    out.query_range.end_ = pattern_pos ;
    out.subject_range.end_ = text_pos ;
    out.transcript.push_terminator();
    out.target_seq = target;
    out.query_source_range = out.query_range;
}
    template<typename It, typename Key>
    Extension::Match target_extension(const Search::Config &cfg, const BlockId id, const Sequence &query, KeyMergeIterator<It,Key> &it){
        std::vector<ExtendedSeed> extended{};
        Extension::Match m = Extension::Match(id, cfg.target->seqs()[id], ::Stats::TargetMatrix(), 0, 0);

        for (auto &hit: it) {
            if ((intersection(hit, extended)))
                continue;

            const Sequence &target = cfg.target->seqs()[id];
            // auto start_trimming = std::chrono::high_resolution_clock::now();
            Sequence query_right = query.subseq(hit.i(), query.length());
            Sequence target_right = target.subseq(hit.j(), target.length());

            std::vector<Letter> query_left = query.subseq(0, hit.i()).reverse();
            std::vector<Letter> target_left = target.subseq(0, hit.j()).reverse();
            //auto end_trimming = std::chrono::high_resolution_clock::now();
            //extend_time.update(3,end_trimming - start_trimming);

            //auto start_ksw = std::chrono::high_resolution_clock::now();
            //Todo: If no cigar is required use ksw2 flag KSW_EZ_EXTZ_ONLY
            KswCigar extension_left = ksw2_align(target_left, query_left, cfg.score_builder->reward(),
                                                 cfg.score_builder->penalty(), cfg.score_builder->gap_open(),
                                                 cfg.score_builder->gap_extend(), 0x40 | 0x80);
           std::cout << " KSW LEFT: ";
           for(auto a: extension_left.cigar_data)
               std::cout << a.first<<a.second;
           std::cout << "\n";

           align_wfa(target_left,query_left,cfg.score_builder->reward(), -cfg.score_builder->penalty(),cfg.score_builder->gap_open(), cfg.score_builder->gap_extend());

            KswCigar extension_right = ksw2_align(target_right, query_right, cfg.score_builder->reward(),
                                                  cfg.score_builder->penalty(), cfg.score_builder->gap_open(),
                                                  cfg.score_builder->gap_extend(), 0x40);

            align_wfa(target_right,query_right,cfg.score_builder->reward(), -cfg.score_builder->penalty(),cfg.score_builder->gap_open(), cfg.score_builder->gap_extend());
            std::cout << " KSW RIGHT: ";
            for(auto a: extension_right.cigar_data)
                std::cout << a.first<<a.second;
            std::cout << "\n";



            KswCigar extension = extension_left + extension_right;
            // auto end_ksw = std::chrono::high_resolution_clock::now();

            //extend_time.update(4,end_ksw - start_ksw);
            //auto start_cigar = std::chrono::high_resolution_clock::now();
            Hsp out = Hsp();
            out.score = extension.score();
            out.bit_score = cfg.score_builder->blast_bit_Score(out.score);
            out.evalue = cfg.score_builder->blast_eValue(out.score, query.length());
            if(out.evalue >= config.max_evalue)
                continue;
            //Todo: if no cigar is required skip the cigar building step
            cigar_to_hsp(extension, cfg.target->seqs()[id], query, hit,out);

            m.hsp.push_back(out);
            extended.emplace_back(hit.i() - extension.max_query(),
                                  out.query_range.end_,
                                  hit.j() - extension.max_target(),
                                  out.subject_range.end_);

            //auto end_cigar = std::chrono::high_resolution_clock::now();
            //extend_time.update(5,end_cigar - start_cigar);
        }
        return m;
    }


std::pair<std::vector<Extension::Match>, Extension::Stats> extend(const Search::Config &cfg,const Sequence &query) {

    //auto start_extend = std::chrono::high_resolution_clock::now();

    std::vector<Extension::Match> matches;
    //auto start_lookup = std::chrono::high_resolution_clock::now();
    auto hits = seed_lookup(query, cfg.target->seqs(), cfg.dna_ref_index.get(), cfg.minimizer_window);
    //auto end_lookup = std::chrono::high_resolution_clock::now();
    //extend_time.update(1,end_lookup-start_lookup);


    //auto start_ungapped = std::chrono::high_resolution_clock::now();
    std::for_each(hits.begin(), hits.end(), [&](SeedMatch &hit) {
        calculate_ungapped_scores(hit, cfg.target->seqs()[hit.id()], query);
    });
    //auto end_ungapped = std::chrono::high_resolution_clock::now();
    //extend_time.update(2,end_ungapped - start_ungapped);


    std::sort(hits.begin(), hits.end(),std::greater<>());
    auto it = merge_keys(hits.begin(), hits.end(), [](const SeedMatch &hit1) { return hit1.id(); });
    while (it.good()) {

        Extension::Match m = target_extension(cfg, it.key(), query,it);
        if (!m.hsp.empty())
            matches.push_back(m);
        ++it;
        //for(auto a: extended_hit_positions)
            //std::cout << a.j_min_extended() << " | " << a.j_max_extended() << " | " << a.i_min_extended() << " | " <<  a.i_max_extended() << std::endl;
    }

        //auto end_extend = std::chrono::high_resolution_clock::now();
        //extend_time.update(0,end_extend - start_extend);
        return {matches, Extension::Stats()};

    }
}






