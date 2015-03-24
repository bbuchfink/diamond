WARN=-Wall -Wno-uninitialized
C=gcc
CXX=g++
CFLAGS=-O3 -DNDEBUG
CXXFLAGS=-O3 -DNDEBUG -Iboost/include $(WARN)

ifneq ($(wildcard boost/.),)
	LIBS=boost/lib/libboost_thread.a boost/lib/libboost_system.a boost/lib/libboost_chrono.a boost/lib/libboost_timer.a boost/lib/libboost_iostreams.a boost/lib/libboost_program_options.a
else
	LIBS=-lboost_thread -lboost_system -lboost_timer -lboost_chrono -lboost_iostreams -lboost_program_options
endif

LIBS+=-lpthread -lz

ifneq (,$(findstring Darwin,$(shell uname)))
else
LIBS+=-lrt
endif

GCC_VER_GTE42 := $(shell echo `gcc -dumpversion | cut -f1-2 -d.` \>= 4.2 | bc )
ifeq ($(GCC_VER_GTE42),1)
CXXFLAGS+=-march=native
endif

OBJS = algo/blast/core/blast_encoding.o algo/blast/core/blast_stat.o algo/blast/core/blast_filter.o algo/blast/core/blast_util.o algo/blast/core/blast_message.o algo/blast/core/ncbi_erf.o algo/blast/core/blast_options.o algo/blast/core/ncbi_math.o algo/blast/core/blast_program.o algo/blast/core/ncbi_std.o algo/blast/core/blast_psi_priv.o algo/blast/core/raw_scoremat.o algo/blast/core/blast_query_info.o algo/blast/core/blast_seg.o main.o basic/options.o util/tinythread.o

DEPS = Makefile align/align_queries.h  align/align_read.h  align/align_sequence.h  align/link_segments.h  align/match_func.h \
basic/const.h		 basic/match.h	     basic/packed_sequence.h	basic/score_matrix.h  basic/sequence.h	basic/shape_config.h  basic/value.h \
basic/edit_transcript.h  basic/options.h     basic/packed_transcript.h	basic/score_traits.h  basic/setup.h	basic/statistics.h    basic/value_type.h \
basic/exceptions.h	 basic/packed_loc.h  basic/reduction.h		basic/seed.h	      basic/shape.h	basic/translate.h \
output/daa_file.h    output/daa_write.h    output/output_buffer.h  output/output_format.h \
output/daa_record.h  output/join_blocks.h  output/output_file.h    output/view.h \
util/async_buffer.h   util/complexity_filter.h	util/fast_mutex.h     util/log_stream.h  util/seq_file_format.h  util/temp_file.h    util/tinythread.h \
util/binary_buffer.h  util/direction.h		util/hash_function.h  util/map.h	 util/system.h		 util/text_buffer.h  util/util.h \
util/binary_file.h    util/double_buffer.h	util/hash_table.h     util/merge_sort.h  util/task_queue.h	 util/thread.h \
data/load_seqs.h  data/queries.h  data/reference.h  data/seed_histogram.h  data/sequence_set.h	data/sorted_list.h  data/string_set.h \
run/make_db.h  run/master_thread.h \
dp/dp_matrix.h	  dp/growing_buffer.h	 dp/scalar_traceback.h	dp/score_vector.h \
dp/floating_sw.h  dp/scalar_dp_matrix.h  dp/score_profile.h	dp/smith_waterman.h \
search/align.h	search/align_range.h  search/align_ungapped.h  search/collision.h  search/hit_filter.h	search/sse_dist.h  search/trace_pt_buffer.h

%.o : %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o : %.c
	$(C) $(CFLAGS) -c -o $@ $<

diamond : $(OBJS)
	$(CXX) -o ../bin/$@ $(OBJS) $(LIBS)

all : diamond
