#include <string.h>
#include "../lib/ksw2/ksw2.h"
#include "ksw2_extension.h"
#include "seed_set_dna.h"


namespace Dna {

// This structure is just for keeping track of extended seeds, also not necessary for the extension
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

// This structure is just for parsing the Ksw Output Cigar to a structure that makes it easier to work with
// This is not necessary for the Wavefront Extension we can also work with the current wavefront output as is
struct KswCigar{
    KswCigar(ksw_extz_t ez):
        max_query_(ez.max_q),
        max_target_(ez.max_t)
        {
        cigar_data.resize(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar_data.emplace_back(ez.cigar[i]>>4, "MID"[ez.cigar[i] & 0xf]);
        }

    KswCigar operator+(const KswCigar& other){
       this->cigar_data.insert(this->cigar_data.end(), other.cigar_data.begin(), other.cigar_data.end());
       return *this;
    }
    int max_query()const{return max_query_;}
    int max_target()const{return max_target_;}
    std::vector<std::pair<int,char>> cigar_data;

private:
    const int max_query_;
    const int max_target_;
};

bool intersection(SeedMatch &hit, const std::vector<ExtendedSeed> &extended){
    return std::any_of(extended.begin(), extended.end(), [hit](ExtendedSeed s)
    {return hit.i() >= s.i_min_extended() && hit.i() + 15 <= s.i_max_extended() && hit.j() >= s.j_min_extended() && hit.j()+15 <= s.j_max_extended() ;});
}

/**
 This is the current ksw2_alignment method.
 Inside we use the "ksw_extz2_sse" method from https://github.com/lh3/ksw2.

 The way we use it, is by trimming the sequences ourselves and extending from the "trimmed" ends of the seeds.
 The parameters for the z-drop in line 79 are not determined for sure  yet, but they give a rough estimate of how the heuristic should behave.

 In the end we are trying to have the same results as BlastN so an x-drop heuristic would of course be desirable to halt the extension.

 */

KswCigar ksw2_align(const Sequence& tseq, const Sequence& qseq, int sc_mch, int sc_mis, int gapo, int gape,int flag) {
    int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = {static_cast<int8_t>(a), static_cast<int8_t>(b), static_cast<int8_t>(b),
                      static_cast<int8_t>(b), 0, static_cast<int8_t>(b), static_cast<int8_t>(a),
                      static_cast<int8_t>(b), static_cast<int8_t>(b), 0, static_cast<int8_t>(b),
                      static_cast<int8_t>(b), static_cast<int8_t>(a), static_cast<int8_t>(b), 0,
                      static_cast<int8_t>(b), static_cast<int8_t>(b), static_cast<int8_t>(b),
                      static_cast<int8_t>(a), 0, 0, 0, 0, 0, 0};


    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));

    ksw_extz2_sse(nullptr, qseq.length(), reinterpret_cast<const uint8_t *>(qseq.data()), tseq.length(),reinterpret_cast<const uint8_t *>(tseq.data()), 5, mat, (int8_t)gapo, (int8_t)gape, -1, 10, 90, flag, &ez);

    KswCigar out(ez);
    free(ez.cigar);

    return out;
}

// this part is also not important for the extension
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

// This part is for converting our intermediate cigar to the diamond output,
// This is also possible with the current output of Wavefront
Hsp cigar_to_hsp(const KswCigar& cigar, const Stats::Blastn_Score &score_builder, const Sequence &target, const Sequence &query, SeedMatch& hit){
    Hsp out = Hsp(true,0);

    int pattern_pos = hit.i() - cigar.max_query() -1;
    int text_pos = hit.j() - cigar.max_target() -1 ;
    int score = 0;
    out.query_range.begin_ = hit.i() - cigar.max_query() -1;
    out.subject_range.begin_ = hit.j() - cigar.max_target() -1;

    for(auto operation: cigar.cigar_data){
        switch(operation.second){
            case 'M':
                for(int j = 0; j < operation.first; ++j){
                    out.push_match(target[text_pos], query[pattern_pos],true);
                    if (target[text_pos] == query[pattern_pos])
                        score += score_builder.reward();
                    else
                        score += score_builder.penalty();
                    pattern_pos++; text_pos++;
                }
                break;
            case 'D':
                out.push_gap(op_deletion,operation.first,target.data()+operation.first+text_pos);
                score -= (score_builder.gap_open() + (score_builder.gap_extend()*(operation.first)));
                text_pos += operation.first;
                break;
            case 'I':
                out.transcript.push_back(op_insertion, (unsigned)operation.first);
                score -= (score_builder.gap_open() + (score_builder.gap_extend()*(operation.first)));
                pattern_pos += operation.first;
                break;
            default:
                break;
        }

    }


    out.score = score;
    out.query_range.end_ = pattern_pos ;
    out.subject_range.end_ = text_pos ;
    out.transcript.push_terminator();
    out.target_seq = target;
    out.query_source_range = out.query_range;

    return out;

}

std::pair<std::vector<Extension::Match>, Extension::Stats> extend(const Search::Config &cfg,const Sequence &query){
    std::vector<Extension::Match> matches;

    auto hits = seed_lookup(query,cfg.target->seqs(),cfg.dna_ref_index.get(),cfg.minimizer_window);
    std::for_each(hits.begin(),hits.end(), [&](SeedMatch &hit){
        calculate_ungapped_scores(hit, cfg.target->seqs()[hit.id()],query);
    });

    std::sort(hits.begin(), hits.end(), [](const SeedMatch &hit1,const SeedMatch &hit2){return hit1.id() < hit2.id();});

    for(BlockId id = 0 ; id < cfg.target->seqs().size(); ++id) {

        auto id_index_start = std::find_if(std::begin(hits), std::end(hits),[&id](SeedMatch const &hit) { return hit.id() == id; });
        if(id_index_start == hits.end())
            continue;
        auto id_index_end = std::find_if(std::begin(hits), std::end(hits),[&id](SeedMatch const &hit) { return hit.id() > id; });

        sort(id_index_start, id_index_end, std::greater<>());
        std::vector<ExtendedSeed> extended_hit_positions{};

        Extension::Match m = Extension::Match(id, cfg.target->seqs()[id], ::Stats::TargetMatrix(), 0,0);

        for (auto hit = id_index_start; hit != id_index_end; ++hit) {
            if ((intersection(*hit, extended_hit_positions)))
                continue;
            const Sequence &target = cfg.target->seqs()[id];
            // here we are simply trimming the sequences to the positions of the seeds
            Sequence query_right = query.subseq(hit->i(), query.length());
            Sequence target_right = target.subseq(hit->j(), target.length());

            std::vector<Letter> query_left = query.subseq(0, hit->i()).reverse();
            std::vector<Letter> target_left = target.subseq(0, hit->j()).reverse();

            KswCigar extension_left = ksw2_align(target_left, query_left, cfg.score_builder->reward(),
                                                 cfg.score_builder->penalty(), cfg.score_builder->gap_open(),
                                                 cfg.score_builder->gap_extend(), 0x40 | 0x80);

            KswCigar extension_right = ksw2_align(target_right, query_right, cfg.score_builder->reward(),
                                                  cfg.score_builder->penalty(), cfg.score_builder->gap_open(),
                                                  cfg.score_builder->gap_extend(), 0x40);


            KswCigar extension = extension_left + extension_right;


            Hsp out = cigar_to_hsp(extension, *cfg.score_builder, cfg.target->seqs()[id], query, *hit);
            out.bit_score = cfg.score_builder->blast_bit_Score(out.score);
            out.evalue = cfg.score_builder->blast_eValue(out.score, query.length());

            if (out.evalue < config.max_evalue) {
                m.hsp.push_back(out);
                extended_hit_positions.emplace_back(hit->i() - extension.max_query(), hit->i() + out.query_range.end_,
                                                    hit->j() - extension.max_target(), hit->j() + out.subject_range.end_);
            }
            if (!m.hsp.empty())
                matches.push_back(m);
        }
    }
    return { matches, Extension::Stats()};

}
}



