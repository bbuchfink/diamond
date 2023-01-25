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
 * VERSION: v20.10.14
 * DESCRIPTION: Simple time profiler
 */

#include "profiler_timer.h"

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

/*
 * System timer
 */
void timer_get_system_time(struct timespec *ts) {
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(),CALENDAR_CLOCK,&cclock);
  clock_get_time(cclock,&mts);
  mach_port_deallocate(mach_task_self(),cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
#else
  clock_gettime(CLOCK_REALTIME,ts);
#endif
}
/*
 * Timers
 */
void timer_start(profiler_timer_t* const timer) {
  timer->accumulated = 0;
  timer_continue(timer);
}
void timer_stop(profiler_timer_t* const timer) {
  timer_pause(timer);
  counter_add(&timer->time_ns,timer->accumulated);
  timer->accumulated = 0;
}
void timer_pause(profiler_timer_t* const timer) {
  struct timespec end_timer;
  timer_get_system_time(&end_timer);
  timer->accumulated += TIME_DIFF_NS(timer->begin_timer,end_timer);
}
void timer_continue(profiler_timer_t* const timer) {
  timer_get_system_time(&timer->begin_timer);
}
void timer_reset(profiler_timer_t* const timer) {
  timer->accumulated = 0;
  counter_reset(&timer->time_ns);
}
uint64_t timer_get_current_lap_ns(profiler_timer_t* const timer) {
  struct timespec end_timer;
  timer_get_system_time(&end_timer);
  return timer->accumulated + TIME_DIFF_NS(timer->begin_timer,end_timer);
}
uint64_t timer_get_current_total_ns(profiler_timer_t* const timer) {
  return counter_get_total(&timer->time_ns) + timer_get_current_lap_ns(timer);
}
uint64_t timer_get_total_ns(const profiler_timer_t* const timer) {
  return counter_get_total(&timer->time_ns) + timer->accumulated;
}
uint64_t timer_get_num_samples(const profiler_timer_t* const timer) {
  return counter_get_num_samples(&timer->time_ns);
}
uint64_t timer_get_min_ns(const profiler_timer_t* const timer) {
  return counter_get_min(&timer->time_ns);
}
uint64_t timer_get_max_ns(const profiler_timer_t* const timer) {
  return counter_get_max(&timer->time_ns);
}
uint64_t timer_get_mean(const profiler_timer_t* const timer) {
  return counter_get_mean(&timer->time_ns);
}
uint64_t timer_get_variance(const profiler_timer_t* const timer) {
  return counter_get_variance(&timer->time_ns);
}
uint64_t timer_get_stddev(const profiler_timer_t* const timer) {
  return counter_get_stddev(&timer->time_ns);
}
void timer_print_total(
    FILE* const stream,
    const profiler_timer_t* const timer) {
  const uint64_t total_time_ns = timer_get_total_ns(timer);
  // Print Total
  if (total_time_ns >= 60000000000ull) {
    fprintf(stream,"%7.2f m ",TIMER_CONVERT_NS_TO_M(total_time_ns));
  } else if (total_time_ns >= 1000000000) {
    fprintf(stream,"%7.2f s ",TIMER_CONVERT_NS_TO_S(total_time_ns));
  } else if (total_time_ns >= 1000000) {
    fprintf(stream,"%7.2f ms",TIMER_CONVERT_NS_TO_MS(total_time_ns));
  } else if (total_time_ns >= 1000) {
    fprintf(stream,"%7.2f us",TIMER_CONVERT_NS_TO_US(total_time_ns));
  } else {
    fprintf(stream,"%7" PRIu64 " ns",total_time_ns);
  }
}
void timer_print(
    FILE* const stream,
    const profiler_timer_t* const timer,
    const profiler_timer_t* const ref_timer) {
  const uint64_t total_time_ns = timer_get_total_ns(timer);
  // Print Total
  timer_print_total(stream,timer);
  // Print percentage wrt reference
  if (ref_timer!=NULL) {
    if (total_time_ns==0) {
        fprintf(stream," (  0.00 %%)");
    } else {
      const uint64_t total_ref_time_ns = timer_get_total_ns(ref_timer);
      if (total_ref_time_ns==0) {
        fprintf(stream," (  n/a  %%)");
      } else {
        const double percentage = (double)(total_time_ns*100)/(double)total_ref_time_ns;
        fprintf(stream," (%6.02f %%)",percentage);
      }
    }
  }
  // Print Calls
  const uint64_t num_calls = timer_get_num_samples(timer);
  if (num_calls >= 1000000000) {
    fprintf(stream," (%5" PRIu64 " Gcalls",num_calls/1000000000);
  } else if (num_calls >= 1000000) {
    fprintf(stream," (%5" PRIu64 " Mcalls",num_calls/1000000);
  } else if (num_calls >= 1000) {
    fprintf(stream," (%5" PRIu64 " Kcalls",num_calls/1000);
  } else if (num_calls > 1 || num_calls == 0) {
    fprintf(stream," (%5" PRIu64 "  calls",num_calls);
  } else {
    fprintf(stream," (%5" PRIu64 "   call",num_calls);
  }
  // Print time/call
  if (num_calls==0) {
    fprintf(stream,",   n/a   s/call)\n");
    return;
  } else {
    const uint64_t ns_per_call = total_time_ns / num_calls;
    if (ns_per_call > 1000000000) {
      fprintf(stream,",%7.2f  s/call",TIMER_CONVERT_NS_TO_S(ns_per_call));
    } else if (ns_per_call > 1000000) {
      fprintf(stream,",%7.2f ms/call",TIMER_CONVERT_NS_TO_MS(ns_per_call));
    } else if (ns_per_call > 1000) {
      fprintf(stream,",%7.2f us/call",TIMER_CONVERT_NS_TO_US(ns_per_call));
    } else {
      fprintf(stream,",%7" PRIu64 " ns/call",ns_per_call);
    }
  }
  // Print Max
  const uint64_t min_ns = timer_get_min_ns(timer);
  if (min_ns > 1000000000) {
    fprintf(stream," {min%.2fs",TIMER_CONVERT_NS_TO_S(min_ns));
  } else if (min_ns > 1000000) {
    fprintf(stream," {min%.2fms",TIMER_CONVERT_NS_TO_MS(min_ns));
  } else if (min_ns > 1000) {
    fprintf(stream," {min%.2fus",TIMER_CONVERT_NS_TO_US(min_ns));
  } else {
    fprintf(stream," {min%" PRIu64 "ns",min_ns);
  }
  // Print Min
  const uint64_t max_ns = timer_get_max_ns(timer);
  if (max_ns > 1000000000) {
    fprintf(stream,",Max%.2fs})\n",TIMER_CONVERT_NS_TO_S(max_ns));
  } else if (max_ns > 1000000) {
    fprintf(stream,",Max%.2fms})\n",TIMER_CONVERT_NS_TO_MS(max_ns));
  } else if (max_ns > 1000) {
    fprintf(stream,",Max%.2fus})\n",TIMER_CONVERT_NS_TO_US(max_ns));
  } else {
    fprintf(stream,",Max%" PRIu64 "ns})\n",max_ns);
  }
}
