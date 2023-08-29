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
 * DESCRIPTION: Wavefront Alignment Algorithms benchmarking tool
 */

#include <omp.h>

#include "align_benchmark_params.h"

#include "utils/commons.h"
#include "utils/sequence_buffer.h"
#include "system/profiler_timer.h"

#include "alignment/score_matrix.h"
#include "edit/edit_dp.h"
#include "gap_linear/nw.h"
#include "gap_affine/swg.h"
#include "gap_affine2p/affine2p_matrix.h"
#include "gap_affine2p/affine2p_dp.h"
#include "wavefront/wavefront_align.h"

#include "benchmark/benchmark_indel.h"
#include "benchmark/benchmark_edit.h"
#include "benchmark/benchmark_gap_linear.h"
#include "benchmark/benchmark_gap_affine.h"
#include "benchmark/benchmark_gap_affine2p.h"

/*
 * WFA lambda (custom match function)
 */
typedef struct {
  char* pattern;
  int pattern_length;
  char* text;
  int text_length;
} match_function_params_t;
match_function_params_t lambda_params;
// Simplest Extend-matching function (for testing purposes)
int lambda_function(int v,int h,void* arguments) {
  // Extract parameters
  match_function_params_t* const match_arguments = (match_function_params_t*)arguments;
  // Check match
  if (v >= match_arguments->pattern_length || h >= match_arguments->text_length) return 0;
  return (match_arguments->pattern[v] == match_arguments->text[h]);
}
/*
 * Algorithms
 */
bool align_benchmark_is_wavefront(
    const alignment_algorithm_type algorithm) {
  return algorithm == alignment_indel_wavefront ||
         algorithm == alignment_edit_wavefront ||
         algorithm == alignment_gap_linear_wavefront ||
         algorithm == alignment_gap_affine_wavefront ||
         algorithm == alignment_gap_affine2p_wavefront;
}
/*
 * Benchmark UTest
 */
void align_pairwise_test() {
  // Patters & Texts
  char * pattern = "GATTACA";
  char * text = "GATCACTA";

  // Penalties
  linear_penalties_t linear_penalties = {
      .match = 0,
      .mismatch = 4,
      .indel = 2,
  };
  affine_penalties_t affine_penalties = {
      .match = 0,
      .mismatch = 4, //9,
      .gap_opening = 6, //13,
      .gap_extension = 2,
  };
  // Ends
  const int pattern_begin_free = 0;
  const int pattern_end_free = 0;
  const int text_begin_free = 0;
  const int text_end_free = 0;
  const bool endsfree =
      pattern_begin_free>0 || pattern_end_free>0 ||
      text_begin_free>0 || text_end_free>0;
  /*
   * Gap-Affine
   */
  // Allocate
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine;
  attributes.linear_penalties = linear_penalties;
  attributes.affine_penalties = affine_penalties;
  attributes.heuristic.strategy = wf_heuristic_none;
  attributes.heuristic.min_wavefront_length = 256;
  attributes.heuristic.max_distance_threshold = 4096;
  attributes.heuristic.steps_between_cutoffs = 10;
  attributes.alignment_scope = compute_alignment; // compute_score
  attributes.memory_mode = wavefront_memory_med;
  attributes.alignment_form.span = (endsfree) ? alignment_endsfree : alignment_end2end;
  attributes.alignment_form.pattern_begin_free = pattern_begin_free;
  attributes.alignment_form.pattern_end_free = pattern_end_free;
  attributes.alignment_form.text_begin_free = text_begin_free;
  attributes.alignment_form.text_end_free = text_end_free;
  attributes.plot.enabled = false;
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  // Align
  wavefront_align(wf_aligner,
      pattern,strlen(pattern),text,strlen(text));
  // CIGAR
  fprintf(stderr,">> WFA2");
  cigar_print_pretty(stderr,wf_aligner->cigar,pattern,strlen(pattern),text,strlen(text));
  fprintf(stderr,"SCORE: %d \n",cigar_score_gap_affine(
      wf_aligner->cigar,&affine_penalties));
  // Plot
  if (attributes.plot.enabled) {
    FILE* const wf_plot = fopen("test.wfa","w");
    wavefront_plot_print(wf_plot,wf_aligner);
    fclose(wf_plot);
  }
  // Free
  wavefront_aligner_delete(wf_aligner);
}
/*
 * Configuration
 */
wavefront_aligner_t* align_input_configure_wavefront(
    align_input_t* const align_input) {
  // Set attributes
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.memory_mode = parameters.wfa_memory_mode;
  if (parameters.wfa_score_only) {
    attributes.alignment_scope = compute_score;
  }
  // WF-Heuristic
  switch (parameters.wfa_heuristic) {
    case wf_heuristic_none:
      attributes.heuristic.strategy = wf_heuristic_none;
      break;
    case wf_heuristic_banded_static:
      attributes.heuristic.strategy = wf_heuristic_banded_static;
      attributes.heuristic.min_k = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_k = parameters.wfa_heuristic_p2;
      break;
    case wf_heuristic_banded_adaptive:
      attributes.heuristic.strategy = wf_heuristic_banded_adaptive;
      attributes.heuristic.min_k = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_k = parameters.wfa_heuristic_p2;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p3;
      break;
    case wf_heuristic_wfadaptive:
      attributes.heuristic.strategy = wf_heuristic_wfadaptive;
      attributes.heuristic.min_wavefront_length = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_distance_threshold = parameters.wfa_heuristic_p2;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p3;
      break;
    case wf_heuristic_xdrop:
      attributes.heuristic.strategy = wf_heuristic_xdrop;
      attributes.heuristic.xdrop = parameters.wfa_heuristic_p1;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p2;
      break;
    case wf_heuristic_zdrop:
      attributes.heuristic.strategy = wf_heuristic_zdrop;
      attributes.heuristic.zdrop = parameters.wfa_heuristic_p1;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p2;
      break;
    default:
      break;
  }
  // Select flavor
  switch (parameters.algorithm) {
    case alignment_indel_wavefront:
      attributes.distance_metric = indel;
      break;
    case alignment_edit_wavefront:
      attributes.distance_metric = edit;
      break;
    case alignment_gap_linear_wavefront:
      attributes.distance_metric = gap_linear;
      attributes.linear_penalties = parameters.linear_penalties;
      break;
    case alignment_gap_affine_wavefront:
      attributes.distance_metric = gap_affine;
      attributes.affine_penalties = parameters.affine_penalties;
      break;
    case alignment_gap_affine2p_wavefront:
      attributes.distance_metric = gap_affine_2p;
      attributes.affine2p_penalties = parameters.affine2p_penalties;
      break;
    default:
      return NULL; // No WF selected
      break;
  }
  // Select alignment form
  attributes.alignment_form.span = (parameters.endsfree) ? alignment_endsfree : alignment_end2end;
  // Misc
  attributes.plot.enabled = (parameters.plot != 0);
  attributes.plot.align_level = (parameters.plot < 0) ? -1 : parameters.plot - 1;
  attributes.system.verbose = parameters.verbose;
  attributes.system.max_memory_abort = parameters.wfa_max_memory;
  attributes.system.max_alignment_score = parameters.wfa_max_score;
  attributes.system.max_num_threads = parameters.wfa_max_threads;
  // Allocate
  return wavefront_aligner_new(&attributes);
}
void align_input_configure_global(
    align_input_t* const align_input) {
  // Clear
  benchmark_align_input_clear(align_input);
  // Penalties
  align_input->linear_penalties = parameters.linear_penalties;
  align_input->affine_penalties = parameters.affine_penalties;
  align_input->affine2p_penalties = parameters.affine2p_penalties;
  // Alignment form
  align_input->ends_free = parameters.endsfree;
  // Output
  align_input->output_file = parameters.output_file;
  align_input->output_full = parameters.output_full;
  // MM
  align_input->mm_allocator = mm_allocator_new(BUFFER_SIZE_1M);
  // WFA
  if (align_benchmark_is_wavefront(parameters.algorithm)) {
    if (parameters.wfa_lambda) {
      align_input->wfa_match_funct = lambda_function;
      align_input->wfa_match_funct_arguments = &lambda_params;
    }
    align_input->wf_aligner = align_input_configure_wavefront(align_input);
  } else {
    align_input->wf_aligner = NULL;
  }
  // PROFILE/STATS
  timer_reset(&align_input->timer);
  // DEBUG
  align_input->debug_flags = 0;
  align_input->debug_flags |= parameters.check_metric;
  if (parameters.check_display) align_input->debug_flags |= ALIGN_DEBUG_DISPLAY_INFO;
  if (parameters.check_correct) align_input->debug_flags |= ALIGN_DEBUG_CHECK_CORRECT;
  if (parameters.check_score) align_input->debug_flags |= ALIGN_DEBUG_CHECK_SCORE;
  if (parameters.check_alignments) align_input->debug_flags |= ALIGN_DEBUG_CHECK_ALIGNMENT;
  align_input->check_linear_penalties = &parameters.linear_penalties;
  align_input->check_affine_penalties = &parameters.affine_penalties;
  align_input->check_affine2p_penalties = &parameters.affine2p_penalties;
  align_input->check_bandwidth = parameters.check_bandwidth;
  align_input->verbose = parameters.verbose;
}
void align_input_configure_local(
    align_input_t* const align_input) {
  // Ends-free configuration
  if (parameters.endsfree) {
    const int plen = align_input->pattern_length;
    const int tlen = align_input->text_length;
    align_input->pattern_begin_free = nominal_prop_u32(plen,parameters.pattern_begin_free);
    align_input->pattern_end_free = nominal_prop_u32(plen,parameters.pattern_end_free);
    align_input->text_begin_free = nominal_prop_u32(tlen,parameters.text_begin_free);
    align_input->text_end_free = nominal_prop_u32(tlen,parameters.text_end_free);
    if (align_benchmark_is_wavefront(parameters.algorithm)) {
      wavefront_aligner_set_alignment_free_ends(align_input->wf_aligner,
          align_input->pattern_begin_free,align_input->pattern_end_free,
          align_input->text_begin_free,align_input->text_end_free);
    }
  }
  // Custom extend-match function
  if (parameters.wfa_lambda) {
    lambda_params.pattern = align_input->pattern;
    lambda_params.pattern_length = align_input->pattern_length;
    lambda_params.text = align_input->text;
    lambda_params.text_length = align_input->text_length;
  }
}
void align_benchmark_free(
    align_input_t* const align_input) {
  if (align_input->wf_aligner) wavefront_aligner_delete(align_input->wf_aligner);
  mm_allocator_delete(align_input->mm_allocator);
}
/*
 * I/O
 */
bool align_benchmark_read_input(
    FILE* input_file,
    char** line1,
    char** line2,
    size_t* line1_allocated,
    size_t* line2_allocated,
    const int seqs_processed,
    align_input_t* const align_input) {
  // Parameters
  int line1_length=0, line2_length=0;
  // Read queries
  line1_length = getline(line1,line1_allocated,input_file);
  if (line1_length==-1) return false;
  line2_length = getline(line2,line2_allocated,input_file);
  if (line1_length==-1) return false;
  // Configure input
  align_input->sequence_id = seqs_processed;
  align_input->pattern = *line1 + 1;
  align_input->pattern_length = line1_length - 2;
  align_input->pattern[align_input->pattern_length] = '\0';
  align_input->text = *line2 + 1;
  align_input->text_length = line2_length - 2;
  if (align_input->text[align_input->text_length] == '\n') {
    align_input->text[align_input->text_length] = '\0';
  }
  return true;
}
/*
 * Display
 */
void align_benchmark_print_progress(
    const int seqs_processed) {
  const uint64_t time_elapsed_alg = timer_get_current_total_ns(&parameters.timer_global);
  const float rate_alg = (float)seqs_processed/(float)TIMER_CONVERT_NS_TO_S(time_elapsed_alg);
  fprintf(stderr,"...processed %d reads (alignment = %2.3f seq/s)\n",seqs_processed,rate_alg);
}
void align_benchmark_print_results(
    align_input_t* const align_input,
    const int seqs_processed,
    const bool print_stats) {
  // Print benchmark results
  fprintf(stderr,"[Benchmark]\n");
  fprintf(stderr,"=> Total.reads            %d\n",seqs_processed);
  fprintf(stderr,"=> Time.Benchmark      ");
  timer_print(stderr,&parameters.timer_global,NULL);
  if (parameters.num_threads == 1) {
    fprintf(stderr,"  => Time.Alignment    ");
    timer_print(stderr,&align_input->timer,&parameters.timer_global);
  } else {
    for (int i=0;i<parameters.num_threads;++i) {
      fprintf(stderr,"  => Time.Alignment.Thread.%0d    ",i);
      timer_print(stderr,&align_input[i].timer,&parameters.timer_global);
    }
  }
  // Print Stats
  const bool checks_enabled =
      parameters.check_display || parameters.check_correct ||
      parameters.check_score || parameters.check_alignments;
  if (checks_enabled && parameters.num_threads==1) {
    const bool print_wf_stats = (parameters.algorithm == alignment_gap_affine_wavefront);
    benchmark_print_stats(stderr,align_input,print_wf_stats);
  }
}
void align_benchmark_plot_wf(
    align_input_t* const align_input,
    const int seq_id) {
  // Setup filename
  char filename[500];
  if (parameters.output_filename != NULL) {
    sprintf(filename,"%s.%03d.plot",parameters.output_filename,seq_id);
  } else {
    sprintf(filename,"%s.%03d.plot",parameters.input_filename,seq_id);
  }
  // Open file
  FILE* const wf_plot = fopen(filename,"w");
  wavefront_plot_print(wf_plot,align_input->wf_aligner);
  fclose(wf_plot);
}
/*
 * Benchmark
 */
void align_benchmark_run_algorithm(
    align_input_t* const align_input) {
  // Sequence-dependent configuration
  align_input_configure_local(align_input);
  // Select algorithm
  switch (parameters.algorithm) {
    // Indel
    case alignment_indel_wavefront:
      benchmark_indel_wavefront(align_input);
      break;
    // Edit
    case alignment_edit_bpm:
      benchmark_edit_bpm(align_input);
      break;
    case alignment_edit_dp:
      benchmark_edit_dp(align_input);
      break;
    case alignment_edit_dp_banded:
      benchmark_edit_dp_banded(align_input,parameters.bandwidth);
      break;
    case alignment_edit_wavefront:
      benchmark_edit_wavefront(align_input);
      break;
    // Gap-linear
    case alignment_gap_linear_nw:
      benchmark_gap_linear_nw(align_input,&parameters.linear_penalties);
      break;
    case alignment_gap_linear_wavefront:
      benchmark_gap_linear_wavefront(align_input,&parameters.linear_penalties);
      break;
    // Gap-affine
    case alignment_gap_affine_swg:
      benchmark_gap_affine_swg(align_input,&parameters.affine_penalties);
      break;
    case alignment_gap_affine_swg_endsfree:
      benchmark_gap_affine_swg_endsfree(
          align_input,&parameters.affine_penalties);
      break;
    case alignment_gap_affine_swg_banded:
      benchmark_gap_affine_swg_banded(align_input,
          &parameters.affine_penalties,parameters.bandwidth);
      break;
    case alignment_gap_affine_wavefront:
      benchmark_gap_affine_wavefront(align_input,&parameters.affine_penalties);
      break;
    // Gap-affine 2p
    case alignment_gap_affine2p_dp:
      benchmark_gap_affine2p_dp(align_input,&parameters.affine2p_penalties);
      break;
    case alignment_gap_affine2p_wavefront:
      benchmark_gap_affine2p_wavefront(align_input,&parameters.affine2p_penalties);
      break;
    default:
      fprintf(stderr,"Algorithm not implemented\n");
      exit(1);
      break;
  }
}
void align_benchmark_sequential() {
  // PROFILE
  timer_reset(&parameters.timer_global);
  timer_start(&parameters.timer_global);
  // I/O files
  parameters.input_file = fopen(parameters.input_filename, "r");
  if (parameters.input_file == NULL) {
    fprintf(stderr,"Input file '%s' couldn't be opened\n",parameters.input_filename);
    exit(1);
  }
  if (parameters.output_filename != NULL) {
    parameters.output_file = fopen(parameters.output_filename, "w");
  }
  // Global configuration
  align_input_t align_input;
  align_input_configure_global(&align_input);
  // Read-align loop
  int seqs_processed = 0, progress = 0;
  while (true) {
    // Read input sequence-pair
    const bool input_read = align_benchmark_read_input(
        parameters.input_file,&parameters.line1,&parameters.line2,
        &parameters.line1_allocated,&parameters.line2_allocated,
        seqs_processed,&align_input);
    if (!input_read) break;
    // Execute the selected algorithm
    align_benchmark_run_algorithm(&align_input);
    // Update progress
    ++seqs_processed;
    if (++progress == parameters.progress) {
      progress = 0;
      if (parameters.verbose >= 0) align_benchmark_print_progress(seqs_processed);
    }
    // DEBUG
    // mm_allocator_print(stderr,align_input.wf_aligner->mm_allocator,false);
    // mm_allocator_print(stderr,align_input.wf_aligner->bialigner->alg_forward->mm_allocator,false);
    // mm_allocator_print(stderr,align_input.wf_aligner->bialigner->alg_reverse->mm_allocator,false);
    // mm_allocator_print(stderr,align_input.wf_aligner->bialigner->alg_subsidiary->mm_allocator,false);
    // Plot
    if (parameters.plot != 0) align_benchmark_plot_wf(&align_input,seqs_processed);
  }
  // Print benchmark results
  timer_stop(&parameters.timer_global);
  if (parameters.verbose >= 0) align_benchmark_print_results(&align_input,seqs_processed,true);
  // Free
  align_benchmark_free(&align_input);
  fclose(parameters.input_file);
  if (parameters.output_file) fclose(parameters.output_file);
  free(parameters.line1);
  free(parameters.line2);
}
void align_benchmark_parallel() {
  // PROFILE
  timer_reset(&parameters.timer_global);
  timer_start(&parameters.timer_global);
  // Open input file
  parameters.input_file = fopen(parameters.input_filename, "r");
  if (parameters.input_file == NULL) {
    fprintf(stderr,"Input file '%s' couldn't be opened\n",parameters.input_filename);
    exit(1);
  }
  if (parameters.output_filename != NULL) {
    parameters.output_file = fopen(parameters.output_filename, "w");
  }
  // Global configuration
  align_input_t align_input[parameters.num_threads];
  for (int tid=0;tid<parameters.num_threads;++tid) {
    align_input_configure_global(align_input+tid);
  }
  // Read-align loop
  sequence_buffer_t* const sequence_buffer = sequence_buffer_new(2*parameters.batch_size,100);
  int seqs_processed = 0, progress = 0, seqs_batch = 0;
  while (true) {
    // Read batch-input sequence-pair
    sequence_buffer_clear(sequence_buffer);
    for (seqs_batch=0;seqs_batch<parameters.batch_size;++seqs_batch) {
      const bool seqs_pending = align_benchmark_read_input(
          parameters.input_file,&parameters.line1,&parameters.line2,
          &parameters.line1_allocated,&parameters.line2_allocated,
          seqs_processed,align_input);
      if (!seqs_pending) break;
      // Add pair pattern-text
      sequence_buffer_add_pair(sequence_buffer,
          align_input->pattern,align_input->pattern_length,
          align_input->text,align_input->text_length);
    }
    if (seqs_batch == 0) break;
    // Parallel processing of the sequences batch
    #pragma omp parallel num_threads(parameters.num_threads)
    {
      int tid = omp_get_thread_num();
      #pragma omp for
      for (int seq_idx=0;seq_idx<seqs_batch;++seq_idx) {
        // Configure sequence
        sequence_offset_t* const offset = sequence_buffer->offsets + seq_idx;
        align_input[tid].sequence_id = seqs_processed;
        align_input[tid].pattern = sequence_buffer->buffer + offset->pattern_offset;
        align_input[tid].pattern_length = offset->pattern_length;
        align_input[tid].text = sequence_buffer->buffer + offset->text_offset;
        align_input[tid].text_length = offset->text_length;
        // Execute the selected algorithm
        align_benchmark_run_algorithm(align_input+tid);
      }
    }
    // Update progress
    seqs_processed += seqs_batch;
    progress += seqs_batch;
    if (progress >= parameters.progress) {
      progress -= parameters.progress;
      if (parameters.verbose >= 0) align_benchmark_print_progress(seqs_processed);
    }
    // DEBUG
    //mm_allocator_print(stderr,align_input->wf_aligner->mm_allocator,false);
  }
  // Print benchmark results
  timer_stop(&parameters.timer_global);
  if (parameters.verbose >= 0) align_benchmark_print_results(align_input,seqs_processed,true);
  // Free
  for (int tid=0;tid<parameters.num_threads;++tid) {
    align_benchmark_free(align_input+tid);
  }
  sequence_buffer_delete(sequence_buffer);
  fclose(parameters.input_file);
  if (parameters.output_file) fclose(parameters.output_file);
  free(parameters.line1);
  free(parameters.line2);
}
/*
 * Main
 */
int main(int argc,char* argv[]) {
  // Parsing command-line options
  parse_arguments(argc,argv);
  // Select option
  if (parameters.algorithm == alignment_test) {
    align_pairwise_test();
  } else {
    // Execute benchmark
    if (parameters.num_threads == 1) {
      align_benchmark_sequential();
    } else {
      align_benchmark_parallel();
    }
  }
}
