#include <iostream>
#include <vector>
#include <cstring>

#include "../sandbox/ksw2.h"

void align(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape)
{
    int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = {
        static_cast<int8_t>(a),static_cast<int8_t>(b),static_cast<int8_t>(b),static_cast<int8_t>(b),0,
        static_cast<int8_t>(b),static_cast<int8_t>(a),static_cast<int8_t>(b),static_cast<int8_t>(b),0,
        static_cast<int8_t>(b),static_cast<int8_t>(b),static_cast<int8_t>(a),static_cast<int8_t>(b),0,
        static_cast<int8_t>(b),static_cast<int8_t>(b),static_cast<int8_t>(b),static_cast<int8_t>(a),0,
        0,0,0,0,0};
    int tl = strlen(tseq), ql = strlen(qseq);
    uint8_t *ts, *qs, c[256];
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];

    ksw_extz2_sse(nullptr, ql, qs, tl, ts, 5, mat, gapo, gape, -1, 30, 100, 0x40 , &ez);
    for (i = 0; i < ez.n_cigar; ++i) {
        printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
    }
    putchar('\n');
    free(ez.cigar); free(ts); free(qs);
}

int main() {
    // This is just a very short example database for the later examples
    const char* target =
        "TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGT"
        "CGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTTGTCCGGGTG"
        "TGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTAC"
        "GTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCC"
        "TCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACTCGAAGGCATT"
        "CAGTACGGTCGTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATACCAGTGGCTTACCGCAAGGTTCTTCTTCGTAAGAACGGTAATA"
        "AAGGAGCTGGTGGCCATAGTTACGGCGCCGATCTAAAGTCATTTGACTTAGGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACTGGAA"
        "CACTAAACATAGCAGTGGTGTTACCCGTGAACTCATGCGTGAGCTTAACGGAGGGGCATACACTCGCTAT";
    int match = 3;
    int mismatch = 3;
    int gapopen = 4;
    int gapextend = 1;

    // Example 1: The query is a perfect subpart of the target, the alignment should end right after the query ends
    const char* query_perfect_match = "TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGC";
    std::cout << "Perfect Match: "; align(target, query_perfect_match, match,mismatch,gapopen, gapextend);

    // Example 2: The query contains 37 Insertions, the alignment should continue over this gap, as there are not enough gaps to trigger the z-drop
    const char* query_insertions = "TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCAAAAAAAAAAACGCGCGCGCGCGCCAAAAAAAAGCGCAGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTTGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAAC";
    std::cout << "Query Insertions: "; align(target, query_insertions, match,mismatch,gapopen, gapextend);

    // Example 3: The query contains a long run of Insertions but only 10 Matches before the insertions, this time the heuristic should be dropped and the alignment should stop after the matches
    const char* query_insertions_long = "TTGTAGATCTAGGGGGGGGCACAGCCTACGCATACATCCCCCCCCCCAAAAAAAAGGGGGGGGGGAAAAAATTTTTTGGGGGGGGAAAAAACCCGCGCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAAC";
    std::cout << "Query long  Insertion: "; align(target, query_insertions_long, match,mismatch,gapopen, gapextend);

    return 0;
}


