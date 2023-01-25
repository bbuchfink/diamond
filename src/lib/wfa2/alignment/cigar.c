/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Edit cigar data-structure (match/mismatch/insertion/deletion)
 */

#include "cigar.h"

/*
 * Setup
 */
cigar_t* cigar_new(
    const int max_operations,
    mm_allocator_t* const mm_allocator) {
  // Allocate
  cigar_t* const cigar = mm_allocator_alloc(mm_allocator,cigar_t);
  // Allocate buffer
  cigar->max_operations = max_operations;
  cigar->operations = mm_allocator_malloc(mm_allocator,cigar->max_operations);
  cigar->begin_offset = 0;
  cigar->end_offset = 0;
  cigar->score = INT32_MIN;
  // MM
  cigar->mm_allocator = mm_allocator;
  // Return
  return cigar;
}
void cigar_clear(
    cigar_t* const cigar) {
  cigar->begin_offset = 0;
  cigar->end_offset = 0;
  cigar->score = INT32_MIN;
}
void cigar_resize(
    cigar_t* const cigar,
    const int max_operations) {
  // Check maximum operations
  if (max_operations > cigar->max_operations) {
    cigar->max_operations = max_operations;
    mm_allocator_free(cigar->mm_allocator,cigar->operations); // Free
    cigar->operations = mm_allocator_malloc(
        cigar->mm_allocator,max_operations); // Allocate
  }
  cigar->begin_offset = 0;
  cigar->end_offset = 0;
  cigar->score = INT32_MIN;
}
void cigar_free(
    cigar_t* const cigar) {
  mm_allocator_free(cigar->mm_allocator,cigar->operations);
  mm_allocator_free(cigar->mm_allocator,cigar);
}
/*
 * Accessors
 */
int cigar_get_matches(
    cigar_t* const cigar) {
  int i, num_matches=0;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    num_matches += (cigar->operations[i]=='M');
  }
  return num_matches;
}
void cigar_add_mismatches(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    cigar_t* const cigar) {
  // Refine adding mismatches
  int i, p=0, t=0;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    // Check limits
    if (p >= pattern_length || t >= text_length) break;
    switch (cigar->operations[i]) {
      case 'M':
        cigar->operations[i] = (pattern[p]==text[t]) ? 'M' : 'X';
        ++p; ++t;
        break;
      case 'I':
        ++t;
        break;
      case 'D':
        ++p;
        break;
      default:
        fprintf(stderr,"[CIGAR] Wrong edit operation\n");
        exit(1);
        break;
    }
  }
  while (p < pattern_length) { cigar->operations[i++] = 'D'; ++p; };
  while (t < text_length) { cigar->operations[i++] = 'I'; ++t; };
  cigar->end_offset = i;
  cigar->operations[cigar->end_offset] = '\0';
  //  // DEBUG
  //  printf("Score=%ld\nPath-length=%" PRIu64 "\nCIGAR=%s\n",
  //      gaba_alignment->score,gaba_alignment->plen,
  //      cigar->operations);
}
/*
 * Score
 */
int cigar_score_edit(
    cigar_t* const cigar) {
  int score = 0, i;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    switch (cigar->operations[i]) {
      case 'M': break;
      case 'X':
      case 'D':
      case 'I': ++score; break;
      default: return INT_MIN;
    }
  }
  return score;
}
int cigar_score_gap_linear(
    cigar_t* const cigar,
    linear_penalties_t* const penalties) {
  int score = 0, i;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    switch (cigar->operations[i]) {
      case 'M': score -= penalties->match; break;
      case 'X': score -= penalties->mismatch; break;
      case 'I': score -= penalties->indel; break;
      case 'D': score -= penalties->indel; break;
      default: return INT_MIN;
    }
  }
  return score;
}
int cigar_score_gap_affine(
    cigar_t* const cigar,
    affine_penalties_t* const penalties) {
  char last_op = '\0';
  int score = 0, i;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    switch (cigar->operations[i]) {
      case 'M':
        score -= penalties->match;
        break;
      case 'X':
        score -= penalties->mismatch;
        break;
      case 'D':
        score -= penalties->gap_extension + ((last_op=='D') ? 0 : penalties->gap_opening);
        break;
      case 'I':
        score -= penalties->gap_extension + ((last_op=='I') ? 0 : penalties->gap_opening);
        break;
      default:
        fprintf(stderr,"[CIGAR] Computing CIGAR score: Unknown operation\n");
        exit(1);
    }
    last_op = cigar->operations[i];
  }
  return score;
}
int cigar_score_gap_affine2p_get_operations_score(
    const char operation,
    const int length,
    affine2p_penalties_t* const penalties) {
  switch (operation) {
    case 'M':
      return penalties->match*length;
    case 'X':
      return penalties->mismatch*length;
    case 'D':
    case 'I': {
      const int score1 = penalties->gap_opening1 + penalties->gap_extension1*length;
      const int score2 = penalties->gap_opening2 + penalties->gap_extension2*length;
      return MIN(score1,score2);
    }
    default:
      fprintf(stderr,"[CIGAR] Computing CIGAR score: Unknown operation\n");
      exit(1);
  }
}
int cigar_score_gap_affine2p(
    cigar_t* const cigar,
    affine2p_penalties_t* const penalties) {
  char last_op = '\0';
  int score = 0, op_length = 0;
  int i;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    // Account for operation
    if (cigar->operations[i] != last_op && last_op != '\0') {
      score -= cigar_score_gap_affine2p_get_operations_score(last_op,op_length,penalties);
      op_length = 0;
    }
    // Add operation
    last_op = cigar->operations[i];
    ++op_length;
  }
  // Account for last operation
  score -= cigar_score_gap_affine2p_get_operations_score(last_op,op_length,penalties);
  return score;
}
/*
 * Utils
 */
int cigar_cmp(
    cigar_t* const cigar_a,
    cigar_t* const cigar_b) {
  // Compare lengths
  const int length_cigar_a = cigar_a->end_offset - cigar_a->begin_offset;
  const int length_cigar_b = cigar_b->end_offset - cigar_b->begin_offset;
  if (length_cigar_a != length_cigar_b) return length_cigar_a - length_cigar_b;
  // Compare operations
  char* const operations_a = cigar_a->operations + cigar_a->begin_offset;
  char* const operations_b = cigar_b->operations + cigar_b->begin_offset;
  int i;
  for (i=0;i<length_cigar_a;++i) {
    if (operations_a[i] != operations_b[i]) {
      return operations_a[i] - operations_b[i];
    }
  }
  // Equal
  return 0;
}
void cigar_copy(
    cigar_t* const cigar_dst,
    cigar_t* const cigar_src) {
  cigar_dst->max_operations = cigar_src->max_operations;
  cigar_dst->begin_offset = cigar_src->begin_offset;
  cigar_dst->end_offset = cigar_src->end_offset;
  cigar_dst->score = cigar_src->score;
  memcpy(cigar_dst->operations+cigar_src->begin_offset,
         cigar_src->operations+cigar_src->begin_offset,
         cigar_src->end_offset-cigar_src->begin_offset);
}
void cigar_append(
    cigar_t* const cigar_dst,
    cigar_t* const cigar_src) {
  // Append
  const int cigar_length = cigar_src->end_offset - cigar_src->begin_offset;
  char* const operations_src = cigar_src->operations + cigar_src->begin_offset;
  char* const operations_dst = cigar_dst->operations + cigar_dst->end_offset;
  memcpy(operations_dst,operations_src,cigar_length);
  // Update offset
  cigar_dst->end_offset += cigar_length;
}
void cigar_append_deletion(
    cigar_t* const cigar,
    const int length) {
  // Append deletions
  char* const operations = cigar->operations + cigar->end_offset;
  int i;
  for (i=0;i<length;++i) {
    operations[i] = 'D';
  }
  // Update offset
  cigar->end_offset += length;
}
void cigar_append_insertion(
    cigar_t* const cigar,
    const int length) {
  // Append insertions
  char* const operations = cigar->operations + cigar->end_offset;
  int i;
  for (i=0;i<length;++i) {
    operations[i] = 'I';
  }
  // Update offset
  cigar->end_offset += length;
}
bool cigar_check_alignment(
    FILE* const stream,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar,
    const bool verbose) {
  // Parameters
  char* const operations = cigar->operations;
  // Traverse CIGAR
  int pattern_pos=0, text_pos=0, i;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    switch (operations[i]) {
      case 'M':
        // Check match
        if (pattern[pattern_pos] != text[text_pos]) {
          if (verbose) {
            fprintf(stream,
                "[AlignCheck] Alignment not matching (pattern[%d]=%c != text[%d]=%c)\n",
                pattern_pos,pattern[pattern_pos],text_pos,text[text_pos]);
          }
          return false;
        }
        ++pattern_pos;
        ++text_pos;
        break;
      case 'X':
        // Check mismatch
        if (pattern[pattern_pos] == text[text_pos]) {
          if (verbose) {
            fprintf(stream,
                "[AlignCheck] Alignment not mismatching (pattern[%d]=%c == text[%d]=%c)\n",
                pattern_pos,pattern[pattern_pos],text_pos,text[text_pos]);
          }
          return false;
        }
        ++pattern_pos;
        ++text_pos;
        break;
      case 'I':
        ++text_pos;
        break;
      case 'D':
        ++pattern_pos;
        break;
      default:
        fprintf(stderr,"[AlignCheck] Unknown edit operation '%c'\n",operations[i]);
        exit(1);
        break;
    }
  }
  // Check alignment length
  if (pattern_pos != pattern_length) {
    if (verbose) {
      fprintf(stream,
          "[AlignCheck] Alignment incorrect length (pattern-aligned=%d,pattern-length=%d)\n",
          pattern_pos,pattern_length);
    }
    return false;
  }
  if (text_pos != text_length) {
    if (verbose) {
      fprintf(stream,
          "[AlignCheck] Alignment incorrect length (text-aligned=%d,text-length=%d)\n",
          text_pos,text_length);
    }
    return false;
  }
  // OK
  return true;
}
/*
 * Display
 */
void cigar_print(
    FILE* const stream,
    cigar_t* const cigar,
    const bool print_matches) {
  // Check null CIGAR
  if (cigar->begin_offset >= cigar->end_offset) return;
  // Print operations
  char last_op = cigar->operations[cigar->begin_offset];
  int last_op_length = 1;
  int i;
  for (i=cigar->begin_offset+1;i<cigar->end_offset;++i) {
    if (cigar->operations[i]==last_op) {
      ++last_op_length;
    } else {
      if (print_matches || last_op != 'M') {
        fprintf(stream,"%d%c",last_op_length,last_op);
      }
      last_op = cigar->operations[i];
      last_op_length = 1;
    }
  }
  if (print_matches || last_op != 'M') {
    fprintf(stream,"%d%c",last_op_length,last_op);
  }
}
int cigar_sprint(
    char* buffer,
    cigar_t* const cigar,
    const bool print_matches) {
  // Parameters
  int pos = 0;
  // Check null CIGAR
  if (cigar->begin_offset >= cigar->end_offset) {
    buffer[pos] = '\0';
    return pos;
  }
  // Print operations
  char last_op = cigar->operations[cigar->begin_offset];
  int last_op_length = 1;
  int i;
  for (i=cigar->begin_offset+1;i<cigar->end_offset;++i) {
    if (cigar->operations[i]==last_op) {
      ++last_op_length;
    } else {
      if (print_matches || last_op != 'M') {
        pos += sprintf(buffer+pos,"%d%c",last_op_length,last_op);
      }
      last_op = cigar->operations[i];
      last_op_length = 1;
    }
  }
  if (print_matches || last_op != 'M') {
    pos += sprintf(buffer+pos,"%d%c",last_op_length,last_op);
  }
  // Return
  buffer[pos] = '\0';
  return pos;
}
void cigar_print_pretty(
    FILE* const stream,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  char* const operations = cigar->operations;
  // Allocate alignment buffers
  const int max_buffer_length = text_length+pattern_length+1;
  char* const pattern_alg = mm_allocator_calloc(mm_allocator,max_buffer_length,char,true);
  char* const ops_alg = mm_allocator_calloc(mm_allocator,max_buffer_length,char,true);
  char* const text_alg = mm_allocator_calloc(mm_allocator,max_buffer_length,char,true);
  // Compute alignment buffers
  int i, alg_pos = 0, pattern_pos = 0, text_pos = 0;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    switch (operations[i]) {
      case 'M':
        if (pattern[pattern_pos] != text[text_pos]) {
          pattern_alg[alg_pos] = pattern[pattern_pos];
          ops_alg[alg_pos] = 'X';
          text_alg[alg_pos++] = text[text_pos];
        } else {
          pattern_alg[alg_pos] = pattern[pattern_pos];
          ops_alg[alg_pos] = '|';
          text_alg[alg_pos++] = text[text_pos];
        }
        pattern_pos++; text_pos++;
        break;
      case 'X':
        if (pattern[pattern_pos] != text[text_pos]) {
          pattern_alg[alg_pos] = pattern[pattern_pos++];
          ops_alg[alg_pos] = ' ';
          text_alg[alg_pos++] = text[text_pos++];
        } else {
          pattern_alg[alg_pos] = pattern[pattern_pos++];
          ops_alg[alg_pos] = 'X';
          text_alg[alg_pos++] = text[text_pos++];
        }
        break;
      case 'I':
        pattern_alg[alg_pos] = '-';
        ops_alg[alg_pos] = ' ';
        text_alg[alg_pos++] = text[text_pos++];
        break;
      case 'D':
        pattern_alg[alg_pos] = pattern[pattern_pos++];
        ops_alg[alg_pos] = ' ';
        text_alg[alg_pos++] = '-';
        break;
      default:
        break;
    }
  }
  i=0;
  while (pattern_pos < pattern_length) {
    pattern_alg[alg_pos+i] = pattern[pattern_pos++];
    ops_alg[alg_pos+i] = '?';
    ++i;
  }
  i=0;
  while (text_pos < text_length) {
    text_alg[alg_pos+i] = text[text_pos++];
    ops_alg[alg_pos+i] = '?';
    ++i;
  }
  // Print alignment pretty
  fprintf(stream,"      ALIGNMENT\t");
  cigar_print(stderr,cigar,true);
  fprintf(stream,"\n");
  fprintf(stream,"      ALIGNMENT.COMPACT\t");
  cigar_print(stderr,cigar,false);
  fprintf(stream,"\n");
  fprintf(stream,"      PATTERN    %s\n",pattern_alg);
  fprintf(stream,"                 %s\n",ops_alg);
  fprintf(stream,"      TEXT       %s\n",text_alg);
  // Free
  mm_allocator_free(mm_allocator,pattern_alg);
  mm_allocator_free(mm_allocator,ops_alg);
  mm_allocator_free(mm_allocator,text_alg);
}


