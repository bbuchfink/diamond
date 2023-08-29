#include <string>
#include "../lib/wfa2/bindings/cpp/WFAligner.hpp"
#include "wfa2_test.h"
#include "../basic/config.h"


namespace WaveExtension {

Hsp cigar_to_hsp(const Stats::Blastn_Score &score_builder, const std::string& cigar,const Sequence &target, const Sequence &query){
    Hsp out = Hsp(true,0);

    int pattern_pos = 0, text_pos = 0, deletion_count = 0, insertion_count = 0, score = 0;
    bool insertion=false, deletion = false;
    out.query_range.begin_ = 0  ;
    out.subject_range.begin_ = 0 ;
    for (auto a: cigar) {
        switch (a) {
            case 'M':
            case 'X':
                if(deletion){
                    out.push_gap(op_deletion,deletion_count,target.data()+deletion_count+text_pos);
                    score -= (score_builder.gap_open() + (score_builder.gap_extend()*(deletion_count-1)));
                    deletion = false;
                    deletion_count = 0;
                }
                if(insertion){
                    out.transcript.push_back(op_insertion, (unsigned)insertion_count);
                    score -= (score_builder.gap_open() + (score_builder.gap_extend()*(insertion_count-1)));
                    insertion = false;
                    insertion_count = 0;
                }
                out.push_match(target[text_pos], query[pattern_pos],true);
                if (target[text_pos] == query[pattern_pos])
                    score += score_builder.reward();
                else
                    score += score_builder.penalty();
                pattern_pos++; text_pos++;
                break;
            case 'D':
                insertion = true;
                insertion_count += 1;
                pattern_pos ++;
                break;
            case 'I':
                deletion = true;
                deletion_count += 1;
                text_pos ++;
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
    out.bit_score = score_builder.blast_bit_Score(out.score);
    out.evalue = score_builder.blast_eValue(out.score,query.length());
    return out;
}

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



void calculate_ungapped_scores(MinimizerHit &hit,const Sequence &target, const Sequence &query){
    int score = 0;
    int i = 0;

    while(query[hit.iMin() - i] == target[hit.jMin() - i]){
        score++;
        i++;
    }

    i = 1 ;
    while(query[hit.iMin() + i] == target[hit.jMin() + i]){
        score++;
        i++;
    }
    hit.score(score);
}

bool intersection(const MinimizerHit &hit,const std::vector<ExtendedSeed> &extended){
    return std::any_of(extended.begin(), extended.end(), [hit](ExtendedSeed s)
        {return hit.iMin() >= s.i_min_extended() && hit.iMax() <= s.i_max_extended() && hit.jMin() >= s.j_min_extended() && hit.jMax() <= s.j_max_extended() ;});
}

std::pair<std::vector<Extension::Match>, Extension::Stats> extend(const Search::Config &cfg,BlockId query_id){
    SequenceSet& target_seqs = cfg.target->seqs();
    std::vector<Extension::Match> matches;
    wfa::WFAlignerGapAffine aligner(-cfg.score_builder->penalty(),cfg.score_builder->gap_open(),cfg.score_builder->gap_extend(),wfa::WFAligner::Alignment,wfa::WFAligner::MemoryHigh);


    for(int i = 0; i < target_seqs.size(); i ++){
        Sequence target_sequence = target_seqs[i];
        std::vector<MinimizerHit> hits = mainMap(cfg,query_id,target_sequence);
        Extension::Match m = Extension::Match(i,target_sequence,::Stats::TargetMatrix(),0,0);

        std::for_each(hits.begin(),hits.end(), [&]( MinimizerHit &hit){
            calculate_ungapped_scores(hit,target_sequence,cfg.query->seqs()[query_id]);
        });

        sort(hits.begin(),hits.end(), std::greater<>());
        std::string query = cfg.query->seqs()[query_id].to_string();
        std::vector<ExtendedSeed> extended_hit_positions;

        for(auto hit: hits){
            if ((intersection(hit,extended_hit_positions)))
                continue;
            std::string target = target_sequence.to_string();
            aligner.alignEndsFree(query,hit.iMin(),hit.iMax(),target,hit.jMin(),hit.jMax()); // Align
            std::string cigar = aligner.getAlignmentCigar();
            Hsp out = cigar_to_hsp(*cfg.score_builder,cigar,target_sequence,cfg.query->seqs()[query_id]);
            if(out.evalue < config.max_evalue){
                m.hsp.push_back(out);
                extended_hit_positions.emplace_back(hit.iMin(),hit.iMax(),hit.jMin(),hit.jMax());
            }

        }
        if(!m.hsp.empty())
            matches.push_back(m);
    }
    return { matches, Extension::Stats()};
}
}