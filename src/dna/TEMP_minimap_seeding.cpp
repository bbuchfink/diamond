#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <map>
#include "../data/block/block.h"
#include "../data/sequence_file.h"
#include "TEMP_minimap_structures.h"
#include "TEMP_minimap_seeding.h"


/**
 *
 * @param char n
 * @return char complement: the complementary nucleotide
 */
char complement(const char& n)
{
    switch(n)
    {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        default:
            break;
    }
    assert(true) ;return (' ');
}

/**
 *  Returns the reverse complement of a DNA String
 * @param std::string sequence
 * @return std::string reverse_complement
 */
std::string reverse_complement(std::string sequence){
    std::reverse(sequence.begin(), sequence.end());

    transform(begin(sequence),end(sequence),begin(sequence),complement);

    return sequence;
}

/**
 *
 * @param std::string s: sequence
 * @param int i : position
 * @param int k : k-mer size
 * @param int r : reverse( if 1 return reverse complement)
 * @return std::string: k-mer of size k at position i
 */
std::string sk(const std::string& s, unsigned long i, int k , int r){

    std::string subsequence = s.substr(i,k);

    if (r)
        return reverse_complement(subsequence);
    else
        return subsequence;
}
/**
 * returns the specific hashValue of a nucleotide
 * @param char: n {A,C,G,T}
 * @return char : {0,1,2,3}
 */
int hash_values(const char& n){
    switch(n)
    {
        case 'A':
            return '0';
        case 'C':
            return '1' ;
        case 'G':
            return '2';
        case 'T':
            return '3';
        default:
            break;

    }
    assert(true) ;return (' ');

}
/**
 *
 * @param std::string s: sequence
 * @return int h: the hashValue of a sequence
 */
int h(const std::string& s){

    std::string  hashedString;

    // transform nucleotide chars to hash Value
    transform(begin(s), end(s), std::back_inserter(hashedString), hash_values);

    //Compute HashValue
    int hashValue = 0;
    for (unsigned long i = 0; i < hashedString.length(); i++)
        hashValue += (hashedString[i] - 48) *pow(4, hashedString.length() - (i + 1));

    return hashValue;
}


/**
 * returns all minimizer of a given sequence s
 * @param std::string s: the DNA sequence from which minimizers are to be extracted
 * @param int k : the k-mer size
 * @param int w : the window-size
 * @return std::vector<minimizer> : dynamic array of minimizer {hashValue, position, strand}
 */
std::unordered_set<Minimizer,MyHashFunction> minimizer_sketch(const std::string& s, const int& w, const int& k ){

    std::unordered_set<Minimizer,MyHashFunction> M;

    for(unsigned long i = 1; i <= s.length()-w-k+1; i++){
        int m = std::numeric_limits<int>::max();

        for( int j = 0; j < w; j++) {
            int u = h(sk(s, i + j, k, 0));
            int v = h(sk(s, i + j, k, 1));

            if (u != v)
                m = std::min({m,u,v},std::less<>());
        }
        for ( int j = 0; j < w; j++){
            int u = h(sk(s, i+j,k,0));
            int v = h(sk(s, i+j,k,1));

            if (u < v && u ==m)
                M.insert(Minimizer(m,i+j,0));
            else if ( v < u && v ==m)
                M.insert(Minimizer(m,i+j,1));
        }
    }


    return M ;
}
/**
 * maps all the minimizers of a given Target-Sequence to their HashValues
 * @param std::vector<fastAentry> fastA : list of Targets
 * @param int w : window-size
 * @param int k : k-mer size
 * @return a HashMap, mapping hashValues to hashIndices {targetSequence, position, Strand}
 */
std::map<int, std::vector<HashIndex>> index(const std::string& target, const int& w, const int& k ){

    std::map<int, std::vector<HashIndex>> H;

    int i = 0;

    std::unordered_set<Minimizer,MyHashFunction> M = minimizer_sketch(target, w,k);
    for( Minimizer m: M){
        if (H.find(m.hashValue) != H.end()){
            if(H[m.hashValue].back().position != m.position ){
                H[m.hashValue].emplace_back(HashIndex(m.position, m.strand));
            }
        }
        else
            H[m.hashValue].emplace_back(HashIndex(m.position, m.strand));
    }


    return H;
}


/**
 *
 * @param map H : a map, mapping hashValues to hashIndeces of the desired target
 * @param std::string q : the given query sequence
 * @param int k : k-mer size
 * @param int w : window size
 * @param int epsilon : the maximum band-width
 * @return vector<mappingResult: a dynamic array containing all successfully mapped queries to the target
 */
std::vector <MinimizerHit> map(const std::map<int, std::vector<HashIndex>>& H, const std::string& q, const int& k, const int& w, const int& epsilon){

    std::vector <MinimizerHit> _result;
    std::vector <MappingMatch> _matches;

    std::unordered_set<Minimizer,MyHashFunction> M = minimizer_sketch(q,w,k);


    for(Minimizer m: M ){
        try{
            const std::vector<HashIndex>& h = H.at(m.hashValue);

            for(HashIndex hI: h){
                if(hI.strand == m.strand) // minimizer on same strand
                    _matches.emplace_back(MappingMatch(0, m.position - hI.position, hI.position));


                else // minimizer on different strand
                    _matches.emplace_back(MappingMatch(1,m.position + hI.position, hI.position));
            }
        }
        catch (const std::out_of_range& oor) {} ;
    }

    sort(_matches.begin(),_matches.end(),comparer());

    unsigned long b = 0;

    for(unsigned long e = 0; e < _matches.size(); e++){

        if(e == _matches.size()-1 || _matches[e+1].strand != _matches[e].strand || _matches[e+1].start - _matches[e].start >= epsilon){

            if(_matches[e].strand == 0)
                _result.emplace_back(MinimizerHit(_matches[b].start + _matches[b].end,
                                                  _matches[b].end,
                                                  (_matches[e].start + _matches[e].end + k) - (_matches[b].start + _matches[b].end)));
            b  = e+1;
        }
    }

    return _result;
};

char const* convert_to_nucl(signed char a){
    switch (a) {
        case 0:
            return "A";
        case 1:
            return "C";
        case 2:
            return "G";
        case 3:
            return "T";
        default:
            return "";
    }
}

std::string sequence_to_string(const Sequence& seq ){
    std::string translated;
    for(int l = 0; l<seq.length();l++){
        translated += convert_to_nucl(seq[l]);
    }
    return translated;
}


std::vector<MinimizerHit> mainMap(const Search::Config &cfg,BlockId query_id,const Sequence &target){
    std::string query = cfg.query->seqs()[query_id].to_string();


    auto H = index(target.to_string(),10,15);

    auto  matches = map(H,query,15,10,1);
/**
 * for(minimizerHit mR: matches){
        std::cout << "match: target= " << mR.targetSequence +1 << " range:"
                  << mR.iMin+1  << "-" << mR.iMax+1 << " to " << mR.jMin+1  << "-"
                  << mR.jMax+1  << " reverse: " << mR.strand << std::endl;


        std::cout << mR.iMax - mR.iMin << " : " << mR.jMax - mR.jMin << std::endl;
    }
 */


    return matches;

}


