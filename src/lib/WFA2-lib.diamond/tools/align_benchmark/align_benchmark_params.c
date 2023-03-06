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
 */

#include "align_benchmark_params.h"

/*
 * Default parameters
 */
align_bench_params_t parameters = {
  // Algorithm
  .algorithm = alignment_edit_wavefront,
  // I/O
  .input_filename = NULL,
  .output_filename = NULL,
  .output_full = false,
  .output_file = NULL,
  // I/O internals
  .input_file = NULL,
  .line1 = NULL,
  .line2 = NULL,
  .line1_allocated = 0,
  .line2_allocated = 0,
  // Penalties
  .linear_penalties = {
      .match = 0,
      .mismatch = 4,
      .indel = 2,
  },
  .affine_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening = 6,
      .gap_extension = 2,
  },
  .affine2p_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening1 = 6,
      .gap_extension1 = 2,
      .gap_opening2 = 24,
      .gap_extension2 = 1,
  },
  // Alignment form
  .endsfree = false,
  .pattern_begin_free = 0.0,
  .text_begin_free = 0.0,
  .pattern_end_free = 0.0,
  .text_end_free = 0.0,
  // Wavefront parameters
  .wfa_score_only = false,
  .wfa_heuristic = wf_heuristic_none,
  .wfa_heuristic_p1 = -1,
  .wfa_heuristic_p2 = -1,
  .wfa_heuristic_p3 = -1,
  .wfa_memory_mode = wavefront_memory_high,
  .wfa_max_memory = UINT64_MAX,
  .wfa_max_score = INT_MAX,
  .wfa_max_threads = 1,
  .wfa_lambda = false,
  // Other algorithms parameters
  .bandwidth = -1,
  // Misc
  .check_bandwidth = -1,
  .check_display = false,
  .check_correct = false,
  .check_score = false,
  .check_alignments = false,
  .check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE,
  .plot = 0,
  // System
  .num_threads = 1,
  .batch_size = 10000,
  .progress = 100000,
  .verbose = 0,
};
/*
 * Menu
 */
void usage() {
  fprintf(stderr,
      "USE: ./align_benchmark -a ALGORITHM -i PATH                             \n"
      "      Options::                                                         \n"
      "        [Algorithm]                                                     \n"
      "          --algorithm|a ALGORITHM                                       \n"
      "            [Indel (Longest Common Subsequence)]                        \n"
      "              indel-wfa                                                 \n"
      "            [Edit (Levenshtein)]                                        \n"
      "              edit-bpm                                                  \n"
      "              edit-dp                                                   \n"
      "              edit-dp-banded                                            \n"
      "              edit-wfa                                                  \n"
      "            [Gap-linear (Needleman-Wunsch)]                             \n"
      "              gap-linear-nw                                             \n"
      "              gap-linear-wfa                                            \n"
      "            [Gap-affine (Smith-Waterman-Gotoh)]                         \n"
      "              gap-affine-swg                                            \n"
      "              gap-affine-swg-banded                                     \n"
      "              gap-affine-wfa                                            \n"
      "            [Gap-affine-2pieces (Concave 2-pieces)]                     \n"
      "              gap-affine2p-dp                                           \n"
      "              gap-affine2p-wfa                                          \n"
      "        [Input & Output]                                                \n"
      "          --input|i PATH                                                \n"
      "          --output|o PATH                                               \n"
      "          --output-full PATH                                            \n"
      "        [Penalties & Span]                                              \n"
      "          --linear-penalties|p M,X,I                                    \n"
      "          --affine-penalties|g M,X,O,E                                  \n"
      "          --affine2p-penalties M,X,O1,E1,O2,E2                          \n"
      "          --ends-free P0,Pf,T0,Tf                                       \n"
      "        [Wavefront parameters]                                          \n"
      "          --wfa-score-only                                              \n"
      "          --wfa-memory-mode 'high'|'med'|'low'|'ultralow'               \n"
      "          --wfa-heuristic STRATEGY                                      \n"
      "          --wfa-heuristic-parameters  P1,P2[,P3]                        \n"
      "            [STRATEGY='banded-static']                                  \n"
      "              P1 = minimum-diagonal-band (e.g., -100)                   \n"
      "              P2 = maximum-diagonal-band (e.g., +100)                   \n"
      "            [STRATEGY='banded-adaptive']                                \n"
      "              P1 = minimum-diagonal-band (e.g., -100)                   \n"
      "              P2 = maximum-diagonal-band (e.g., +100)                   \n"
      "              P3 = steps-between-cutoffs                                \n"
      "            [STRATEGY='wfa-adaptive']                                   \n"
      "              P1 = minimum-wavefront-length                             \n"
      "              P2 = maximum-difference-distance                          \n"
      "              P3 = steps-between-cutoffs                                \n"
      "            [STRATEGY='xdrop']                                          \n"
      "              P1 = x-drop                                               \n"
      "              P2 = steps-between-cutoffs                                \n"
      "            [STRATEGY='zdrop']                                          \n"
      "              P1 = z-drop                                               \n"
      "              P2 = steps-between-cutoffs                                \n"
      "          --wfa-max-memory BYTES                                        \n"
      "          --wfa-max-score INT                                           \n"
      "          --wfa-max-threads INT (intra-parallelism; default=1)          \n"
      "        [Other Parameters]                                              \n"
      "          --bandwidth INT                                               \n"
      "        [Misc]                                                          \n"
      "          --check|c 'correct'|'score'|'alignment'                       \n"
      "          --check-distance 'indel'|'edit'|'linear'|'affine'|'affine2p'  \n"
      "          --check-bandwidth INT                                         \n"
      "          --plot                                                        \n"
      "        [System]                                                        \n"
      "          --num-threads|t INT                                           \n"
      "          --batch-size INT                                              \n"
    //"          --progress|P INT                                              \n"
      "          --verbose|v INT                                               \n"
      "          --quiet|q                                                     \n"
      "          --help|h                                                      \n");
}
/*
 * Parsing arguments
 */
void parse_arguments(
    int argc,
    char** argv) {
  struct option long_options[] = {
    /* Input */
    { "algorithm", required_argument, 0, 'a' },
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "output-full", required_argument, 0, 800 },
    /* Penalties */
    { "linear-penalties", required_argument, 0, 'p' },
    { "affine-penalties", required_argument, 0, 'g' },
    { "affine2p-penalties", required_argument, 0, 900 },
    { "ends-free", required_argument, 0, 901 },
    /* Wavefront parameters */
    { "wfa-score-only", no_argument, 0, 1000 },
    { "wfa-memory-mode", required_argument, 0, 1001 },
    { "wfa-heuristic", required_argument, 0, 1002 },
    { "wfa-heuristic-parameters", required_argument, 0, 1003 },
    { "wfa-max-memory", required_argument, 0, 1005 },
    { "wfa-max-score", required_argument, 0, 1006 },
    { "wfa-max-threads", required_argument, 0, 1007 },
    { "wfa-lambda", no_argument, 0, 1008 },
    /* Other alignment parameters */
    { "bandwidth", required_argument, 0, 2000 },
    /* Misc */
    { "check", required_argument, 0, 'c' },
    { "check-distance", required_argument, 0, 3001 },
    { "check-bandwidth", required_argument, 0, 3002 },
    { "plot", optional_argument, 0, 3003 },
    /* System */
    { "num-threads", required_argument, 0, 't' },
    { "batch-size", required_argument, 0, 4000 },
    { "progress", required_argument, 0, 'P' },
    { "verbose", required_argument, 0, 4001 },
    { "verbose1", no_argument, 0, 'v' },
    { "quiet", no_argument, 0, 'q' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  if (argc <= 1) {
    usage();
    exit(0);
  }
  while (1) {
    c=getopt_long(argc,argv,"a:i:o:p:g:P:c:vqt:h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /*
     * Algorithm
     */
    case 'a': {
      /* Test bench */
      if (strcmp(optarg,"test")==0) {
        parameters.algorithm = alignment_test;
      // Indel
      } else if (strcmp(optarg,"indel-wfa")==0) {
        parameters.algorithm = alignment_indel_wavefront;
      // Edit
      } else if (strcmp(optarg,"edit-bpm")==0) {
        parameters.algorithm = alignment_edit_bpm;
      } else if (strcmp(optarg,"edit-dp")==0) {
        parameters.algorithm = alignment_edit_dp;
      } else if (strcmp(optarg,"edit-dp-banded")==0) {
        parameters.algorithm = alignment_edit_dp_banded;
      } else if (strcmp(optarg,"edit-wfa")==0) {
        parameters.algorithm = alignment_edit_wavefront;
      // Gap-Linear
      } else if (strcmp(optarg,"gap-linear-nw")==0 ||
                 strcmp(optarg,"gap-linear-dp")==0) {
        parameters.algorithm = alignment_gap_linear_nw;
      } else if (strcmp(optarg,"gap-linear-wfa")==0) {
        parameters.algorithm = alignment_gap_linear_wavefront;
      // Gap-Affine
      } else if (strcmp(optarg,"gap-affine-swg")==0 ||
                 strcmp(optarg,"gap-affine-dp")==0) {
        parameters.algorithm = alignment_gap_affine_swg;
      } else if (strcmp(optarg,"gap-affine-swg-banded")==0 ||
                 strcmp(optarg,"gap-affine-dp-banded")==0) {
        parameters.algorithm = alignment_gap_affine_swg_banded;
      } else if (strcmp(optarg,"gap-affine-wfa")==0) {
        parameters.algorithm = alignment_gap_affine_wavefront;
      // Gap-Affine 2-Pieces
      } else if (strcmp(optarg,"gap-affine2p-dp")==0) {
        parameters.algorithm = alignment_gap_affine2p_dp;
      } else if (strcmp(optarg,"gap-affine2p-wfa")==0) {
        parameters.algorithm = alignment_gap_affine2p_wavefront;
      } else {
        fprintf(stderr,"Algorithm '%s' not recognized\n",optarg);
        exit(1);
      }
      break;
    }
    /*
     * Input/Output
     */
    case 'i':
      parameters.input_filename = optarg;
      break;
    case 'o':
      parameters.output_filename = optarg;
      break;
    case 800: // --output-full
      parameters.output_filename = optarg;
      parameters.output_full = true;
      break;
    /*
     * Penalties
     */
    case 'p': { // --linear-penalties M,X,I
      char* sentinel = strtok(optarg,",");
      parameters.linear_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.linear_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.linear_penalties.indel = atoi(sentinel);
      break;
    }
    case 'g': { // --affine-penalties M,X,O,E
      char* sentinel = strtok(optarg,",");
      parameters.affine_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.gap_opening = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.gap_extension = atoi(sentinel);
      break;
    }
    case 900: { // --affine2p-penalties M,X,O1,E1,O2,E2
      char* sentinel = strtok(optarg,",");
      parameters.affine2p_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_opening1 = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_extension1 = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_opening2 = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_extension2 = atoi(sentinel);
      break;
    }
    case 901: { // --ends-free P0,Pf,T0,Tf
      parameters.endsfree = true;
      char* sentinel = strtok(optarg,",");
      parameters.pattern_begin_free = atof(sentinel);
      sentinel = strtok(NULL,",");
      parameters.pattern_end_free = atof(sentinel);
      sentinel = strtok(NULL,",");
      parameters.text_begin_free = atof(sentinel);
      sentinel = strtok(NULL,",");
      parameters.text_end_free = atof(sentinel);
      break;
    }
    /*
     * Wavefront parameters
     */
    case 1000: // --wfa-score-only
      parameters.wfa_score_only = true;
      break;
    case 1001: // --wfa-memory-mode in {'high','med','low'}
      if (strcmp(optarg,"high")==0) {
        parameters.wfa_memory_mode = wavefront_memory_high;
      } else if (strcmp(optarg,"med")==0) {
        parameters.wfa_memory_mode = wavefront_memory_med;
      } else if (strcmp(optarg,"low")==0) {
        parameters.wfa_memory_mode = wavefront_memory_low;
      } else if (strcmp(optarg,"ultralow")==0) {
        parameters.wfa_memory_mode = wavefront_memory_ultralow;
      } else {
        fprintf(stderr,"Option '--wfa-memory-mode' must be in {'high','med','low','ultralow'}\n");
        exit(1);
      }
      break;
    case 1002: // --wfa-heuristic in {'none'|'banded-static'|'banded-adaptive'|'wfa-adaptive'|'xdrop'|'zdrop'}
      if (strcmp(optarg,"none")==0) {
        parameters.wfa_heuristic = wf_heuristic_none;
      } else if (strcmp(optarg,"banded-static")==0 || strcmp(optarg,"banded")==0) {
        parameters.wfa_heuristic = wf_heuristic_banded_static;
      } else if (strcmp(optarg,"banded-adaptive")==0) {
        parameters.wfa_heuristic = wf_heuristic_banded_adaptive;
      } else if (strcmp(optarg,"wfa-adaptive")==0) {
        parameters.wfa_heuristic = wf_heuristic_wfadaptive;
      } else if (strcmp(optarg,"xdrop")==0) {
        parameters.wfa_heuristic = wf_heuristic_xdrop;
      } else if (strcmp(optarg,"zdrop")==0) {
        parameters.wfa_heuristic = wf_heuristic_zdrop;
      } else {
        fprintf(stderr,"Option '--wf-heuristic' must be in {'none'|'banded-static'|'banded-adaptive'|'wfa-adaptive'|'xdrop'|'zdrop'}\n");
        exit(1);
      }
      break;
    case 1003: { // --wfa-heuristic-parameters  <P1>,<P2>[,<P3>]
      char* sentinel = strtok(optarg,",");
      const int p1 = atoi(sentinel);
      parameters.wfa_heuristic_p1 = p1;
      sentinel = strtok(NULL,",");
      const int p2 = atoi(sentinel);
      parameters.wfa_heuristic_p2 = p2;
      sentinel = strtok(NULL,",");
      if (sentinel != NULL) {
        const int p3 = atoi(sentinel);
        parameters.wfa_heuristic_p3 = p3;
      }
      break;
    }
    case 1005: // --wfa-max-memory
      parameters.wfa_max_memory = atol(optarg);
      break;
    case 1006: // --wfa-max-score
      parameters.wfa_max_score = atoi(optarg);
      break;
    case 1007: // --wfa-max-threads
      parameters.wfa_max_threads = atoi(optarg);
      break;
    case 1008: // --wfa-lambda
      parameters.wfa_lambda = true;
      break;
    /*
     * Other alignment parameters
     */
    case 2000: // --bandwidth
      parameters.bandwidth = atoi(optarg);
      break;
    /*
     * Misc
     */
    case 'c':
      if (optarg ==  NULL) { // default = correct
        parameters.check_correct = true;
        parameters.check_score = false;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"display")==0) {
        parameters.check_display = true;
      } else if (strcasecmp(optarg,"correct")==0) {
        parameters.check_correct = true;
        parameters.check_score = false;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"score")==0) {
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"alignment")==0) {
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = true;
      } else {
        fprintf(stderr,"Option '--check' must be in {'correct','score','alignment'}\n");
        exit(1);
      }
      break;
    case 3001: // --check-distance in {'indel','edit','linear','affine','affine2p'}
      if (strcasecmp(optarg,"indel")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_INDEL;
      } else if (strcasecmp(optarg,"edit")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_EDIT;
      } else if (strcasecmp(optarg,"linear")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_LINEAR;
      } else if (strcasecmp(optarg,"affine")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE;
      } else if (strcasecmp(optarg,"affine2p")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE2P;
      } else {
        fprintf(stderr,"Option '--check-distance' must be in {'indel','edit','linear','affine','affine2p'}\n");
        exit(1);
      }
      break;
    case 3002: // --check-bandwidth
      parameters.check_bandwidth = atoi(optarg);
      break;
    case 3003: // --plot
      parameters.plot = (optarg==NULL) ? 1 : atoi(optarg);
      break;
    /*
     * System
     */
    case 't': // --num-threads
      parameters.num_threads = atoi(optarg);
      break;
    case 4000: // --batch-size
      parameters.batch_size = atoi(optarg);
      break;
    case 'P':
      parameters.progress = atoi(optarg);
      break;
    case 'v':
      parameters.verbose = 1;
      break;
    case 4001: // --verbose (long option)
      parameters.verbose = atoi(optarg);
      if (parameters.verbose < 0 || parameters.verbose > 4) {
        fprintf(stderr,"Option '--verbose' must be in {0,1,2,3,4}\n");
        exit(1);
      }
      break;
    case 'q':
      parameters.verbose = -1;
      break;
    case 'h':
      usage();
      exit(1);
    // Other
    case '?': default:
      fprintf(stderr,"Option not recognized \n");
      exit(1);
    }
  }
  // Checks general
  if (parameters.algorithm!=alignment_test && parameters.input_filename==NULL) {
    fprintf(stderr,"Option --input is required \n");
    exit(1);
  }
  // Check 'ends-free' parameter
  if (parameters.endsfree) {
    switch (parameters.algorithm) {
      case alignment_gap_affine_swg:
        parameters.algorithm = alignment_gap_affine_swg_endsfree;
        break;
      case alignment_indel_wavefront:
      case alignment_edit_wavefront:
      case alignment_gap_linear_wavefront:
      case alignment_gap_affine_wavefront:
      case alignment_gap_affine2p_wavefront:
        break;
      default:
        fprintf(stderr,"Ends-free variant not implemented for the selected algorithm\n");
        exit(1);
        break;
    }
  }
  // Check 'bandwidth' parameter
  switch (parameters.algorithm) {
    case alignment_edit_dp_banded:
    case alignment_gap_affine_swg_banded:
      if (parameters.bandwidth == -1) {
        fprintf(stderr,"Parameter 'bandwidth' has to be provided for banded algorithms\n");
        exit(1);
      }
      break;
    default:
      if (parameters.bandwidth != -1) {
        fprintf(stderr,"Parameter 'bandwidth' has no effect with the selected algorithm\n");
        exit(1);
      }
      break;
  }
  // Check 'wfa-heuristic'
  switch (parameters.wfa_heuristic) {
    case wf_heuristic_banded_static:
    case wf_heuristic_xdrop:
    case wf_heuristic_zdrop:
      if (parameters.wfa_heuristic_p1 == -1 ||
          parameters.wfa_heuristic_p2 == -1) {
        fprintf(stderr,"Heuristic requires parameters '--wfa-heuristic-parameters' <P1>,<P2>\n");
        exit(1);
      }
      break;
    case wf_heuristic_banded_adaptive:
    case wf_heuristic_wfadaptive:
      if (parameters.wfa_heuristic_p1 == -1 ||
          parameters.wfa_heuristic_p2 == -1 ||
          parameters.wfa_heuristic_p3 == -1) {
        fprintf(stderr,"Heuristic requires parameters '--wfa-heuristic-parameters' <P1>,<P2>,<P3>\n");
        exit(1);
      }
      break;
    default:
      break;
  }
  // Checks parallel
  if (parameters.num_threads > 1) {
    if (parameters.plot > 0) {
      fprintf(stderr,"Parameter 'plot' disabled for parallel executions\n");
      parameters.plot = 0;
    }
  }
}
