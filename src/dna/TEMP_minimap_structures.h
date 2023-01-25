#include <iostream>
#include "../data/block/block.h"
#include "../data/sequence_file.h"

#ifndef STRUCTURES
#define STRUCTURES


struct FastAentry{
    std::string header;
    std::string sequence;
    FastAentry(const std::string& header_ , const std::string& sequence_ ){
        header = header_;
        sequence = sequence_;
    };
};


struct Minimizer{
    int hashValue;
    unsigned long position;
    int strand;

    Minimizer(const int& hashvalue_, const unsigned long& position_, const int& strand_ ){
        hashValue = hashvalue_;
        position = position_;
        strand = strand_ ;
    }
    bool operator==(const Minimizer& t) const
    {
        return (this->hashValue == t.hashValue);
    }

};

class MyHashFunction {
public:
    size_t operator()(const Minimizer& t) const
    {
        return t.hashValue;
    }
};

struct HashIndex{
    int position;
    int strand;

    HashIndex(const int position_ , const int strand_ ){
        position = position_;
        strand = strand_;
    }
};

struct MinimizerHit{


    MinimizerHit(const int i, const int j, const int length):
    i_(i),
    j_(j),
    length_(length)

    {};

    int iMin()const{return i_;}
    int iMax()const{return i_ + length_;}
    int jMin()const{return j_;}
    int jMax()const{return j_ + length_;}
    int length()const{return length_;}
    int score()const{return score_;}
    void score(int score){score_ = score;}
    bool operator>(const MinimizerHit& hit)const
    {return  this->score() > hit.score();}

private:
    int i_,j_,length_, score_;

};


struct MappingMatch{
    int strand;
    int start;
    int end;

    MappingMatch(const int& strand_, const int& start_, const int& end_){
        strand = strand_;
        start = start_;
        end = end_;
    }
};

// structure for sorting the mappingMatches
struct comparer{
    inline bool operator()(const MappingMatch& one, const MappingMatch& two){


        return ((one.strand < two.strand) ||
                (one.strand == two.strand && one.start < two.start) ||
                (one.strand == two.strand && one.start == two.start && one.end < two.end));

    }
};





#endif
