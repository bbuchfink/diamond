cmake_minimum_required (VERSION 2.6)
project (DIAMOND)

option(BUILD_STATIC "BUILD_STATIC" OFF)

if(BUILD_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBRARIES OFF)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCHNATIVE)
if(COMPILER_SUPPORTS_MARCHNATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-uninitialized -Wno-deprecated-declarations")

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
  src/dp/floating_sw.cpp
  src/align/align_sequence_anchored.cpp
  src/align/align_sequence_simple.cpp
  src/basic/hssp.cpp
  src/dp/ungapped_align.cpp
  src/run/tools.cpp
  src/dp/greedy_align.cpp
  src/run/benchmark.cpp
  src/search/stage2.cpp
  src/output/output_format.cpp
)

target_link_libraries(diamond ${ZLIB_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS diamond DESTINATION bin)
