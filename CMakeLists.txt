cmake_minimum_required (VERSION 2.6)
project (DIAMOND)

option(BUILD_STATIC "BUILD_STATIC" OFF)

if(BUILD_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBRARIES OFF)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

if(CMAKE_BUILD_MARCH)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${CMAKE_BUILD_MARCH}")
else()
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCHNATIVE)
  if(COMPILER_SUPPORTS_MARCHNATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif()
endif()

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-uninitialized -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-unused-variable")

include_directories(
  "${CMAKE_SOURCE_DIR}/src"
  "${ZLIB_INCLUDE_DIR}")

add_executable(diamond src/run/main.cpp
  src/basic/config.cpp
  src/util/tinythread.cpp
  src/util/compressed_stream.cpp
  src/basic/score_matrix.cpp
  src/blast/blast_filter.cpp
  src/blast/blast_seg.cpp
  src/blast/sm_blosum45.c
  src/blast/sm_blosum50.c
  src/blast/sm_blosum62.c
  src/blast/sm_blosum80.c
  src/blast/sm_blosum90.c
  src/blast/sm_pam30.c
  src/blast/sm_pam70.c
  src/blast/sm_pam250.c
  src/data/queries.cpp
  src/data/reference.cpp
  src/data/seed_histogram.cpp
  src/output/daa_record.cpp
  src/search/search.cpp
  src/util/command_line_parser.cpp
  src/util/seq_file_format.cpp
  src/util/util.cpp 
  src/util/Timer.cpp
  src/basic/basic.cpp
  src/basic/hssp.cpp
  src/dp/ungapped_align.cpp
  src/run/tools.cpp
  src/dp/greedy_align.cpp
  src/run/benchmark.cpp
  src/search/stage2.cpp
  src/output/output_format.cpp
  src/output/join_blocks.cpp
  src/data/frequent_seeds.cpp
  src/align/query_mapper.cpp
  src/align/align_target.cpp
  src/output/blast_tab_format.cpp
  src/dp/padded_banded_sw.cpp
  src/dp/needleman_wunsch.cpp
  src/output/blast_pairwise_format.cpp
  src/extra/roc.cpp
  src/dp/comp_based_stats.cpp
  src/extra/model_sim.cpp
  src/run/double_indexed.cpp
  src/search/collision.cpp
  src/output/sam_format.cpp
  src/align/align.cpp
  src/search/setup.cpp
  src/extra/opt.cpp
  src/dp/diag_scores.cpp
  src/data/taxonomy.cpp
  src/lib/tantan/tantan.cc
  src/basic/masking.cpp
  src/dp/swipe.cpp
  src/dp/banded_sw.cpp
  src/data/sorted_list.cpp
  src/data/seed_set.cpp
  src/util/binary_file.cpp
  src/util/simd.cpp
  src/output/taxon_format.cpp
  src/output/view.cpp
  src/output/output_sink.cpp
  src/output/target_culling.cpp
)

target_link_libraries(diamond ${ZLIB_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS diamond DESTINATION bin)
