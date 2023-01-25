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
 * VERSION: v20.08.25
 * DESCRIPTION: Simple profile counter
 */

#ifndef PROFILER_COUNTER_H_
#define PROFILER_COUNTER_H_

#include "../utils/commons.h"

/*
 * Counters
 */
typedef struct {
  uint64_t total;
  uint64_t samples;
  uint64_t min;
  uint64_t max;
  double m_oldM;
  double m_newM;
  double m_oldS;
  double m_newS;
} profiler_counter_t;

void counter_reset(
    profiler_counter_t* const counter);
void counter_add(
    profiler_counter_t* const counter,
    const uint64_t amount);

uint64_t counter_get_total(const profiler_counter_t* const counter);
uint64_t counter_get_num_samples(const profiler_counter_t* const counter);
uint64_t counter_get_min(const profiler_counter_t* const counter);
uint64_t counter_get_max(const profiler_counter_t* const counter);
double counter_get_mean(const profiler_counter_t* const counter);
double counter_get_variance(const profiler_counter_t* const counter);
double counter_get_stddev(const profiler_counter_t* const counter);

void counter_combine_sum(
    profiler_counter_t* const counter_dst,
    profiler_counter_t* const counter_src);

void counter_print(
    FILE* const stream,
    const profiler_counter_t* const counter,
    const profiler_counter_t* const ref_counter,
    const char* const units,
    const bool full_report);
void percentage_print(
    FILE* const stream,
    const profiler_counter_t* const counter,
    const char* const units);

/*
 * Reference Counter (Counts wrt a reference counter. Eg ranks)
 */
typedef struct {
  uint64_t begin_count;       // Counter
  profiler_counter_t counter; // Total count & samples taken
  uint64_t accumulated;       // Total accumulated
} profiler_rcounter_t;

void rcounter_start(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference);
void rcounter_stop(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference);
void rcounter_pause(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference);
void rcounter_restart(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference);
void rcounter_reset(
    profiler_rcounter_t* const rcounter);

uint64_t rcounter_get_total(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_num_samples(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_min(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_max(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_mean(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_variance(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_stddev(profiler_rcounter_t* const rcounter);

/*
 * Display
 */
#define PRIcounter "lu(#%" PRIu64 ",m%" PRIu64 ",M%" PRIu64",{%.2f})"
#define PRIcounterVal(counter) \
  counter_get_total(counter), \
  counter_get_num_samples(counter), \
  counter_get_min(counter), \
  counter_get_max(counter), \
  counter_get_mean(counter)
#define PRIcounterX "lu(#%" PRIu64 ",m%" PRIu64 ",M%" PRIu64 ",{%.2f,%.2f,%.2f})"
#define PRIcounterXVal(counter) \
  counter_get_total(counter), \
  counter_get_num_samples(counter), \
  counter_get_min(counter), \
  counter_get_max(counter), \
  counter_get_mean(counter), \
  counter_get_variance(counter), \
  counter_get_stddev(counter)

#endif /* PROFILER_COUNTER_H_ */
