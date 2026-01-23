# Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 2.6...3.5)
project (DIAMOND)
include(CheckCXXCompilerFlag)
include(CheckSymbolExists)
include(CheckTypeSize)

option(BUILD_STATIC "BUILD_STATIC" OFF)
option(STATIC_LIBGCC "STATIC_LIBGCC" OFF)
option(STATIC_LIBSTDC++ "STATIC_LIBSTDC++" OFF)
option(X86 "X86" ON)
option(ARM "ARM" OFF)
option(AARCH64 "AARCH64" OFF)
option(WITH_ZSTD "WITH_ZSTD" OFF)
option(CROSS_COMPILE "CROSS_COMPILE" OFF)

# experimental options not suitable for deployment purposes
option(EXTRA "EXTRA" OFF)
option(STRICT_BAND "STRICT_BAND" ON)
option(LEFTMOST_SEED_FILTER "LEFTMOST_SEED_FILTER" ON)
option(SEQ_MASK "SEQ_MASK" ON)
option(DP_STAT "DP_STAT" OFF)
option(HIT_KEEP_TARGET_ID "HIT_KEEP_TARGET_ID" OFF)
option(LONG_SEEDS "LONG_SEEDS" OFF)
option(WITH_AVX512 "WITH_AVX512" OFF)
option(WITH_DNA "WITH_DNA" OFF)
option(WITH_MCL "WITH_MCL" OFF)
option(USE_TLS "USE_TLS" OFF)

set(MAX_SHAPE_LEN 19)

if(NOT CROSS_COMPILE AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*|arm64.*)")
  set(X86 OFF)
  set(AARCH64 ON)
elseif(NOT CROSS_COMPILE AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm.*|ARM.*)")
  set(X86 OFF)
  set(ARM ON)
elseif(NOT CROSS_COMPILE AND CMAKE_SYSTEM_PROCESSOR MATCHES "PPC64*|ppc64*|powerpc64*")
  set(X86 OFF)
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -maltivec")
elseif(NOT CROSS_COMPILE AND CMAKE_SYSTEM_PROCESSOR MATCHES "^s390|sparc")
  set(X86 OFF)
endif()

if(STATIC_LIBSTDC++)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")
endif()

if(STRICT_BAND)
  add_definitions(-DSTRICT_BAND)
endif()

if(USE_TLS)
  add_definitions(-DUSE_TLS)
endif()

if(SEQ_MASK)
  add_definitions(-DSEQ_MASK)
endif()

if(LEFTMOST_SEED_FILTER)
  add_definitions(-DLEFTMOST_SEED_FILTER)
endif()

if(DP_STAT)
  add_definitions(-DDP_STAT)
endif()

if(EXTRA)
  add_definitions(-DEXTRA)
endif()

if(HIT_KEEP_TARGET_ID)
  add_definitions(-DHIT_KEEP_TARGET_ID)
endif()

if(LONG_SEEDS)
  add_definitions(-DLONG_SEEDS)
endif()

if(WITH_AVX512)
  add_definitions(-DWITH_AVX512)
endif()

if(WITH_DNA)
  add_definitions(-DWITH_DNA)
endif()

if(WITH_MCL)
  add_definitions(-DWITH_MCL)
endif()

if(WITH_FAMSA)
  add_definitions(-DWITH_FAMSA)
endif()

add_definitions(-DMAX_SHAPE_LEN=${MAX_SHAPE_LEN})
add_definitions(-D_ITERATOR_DEBUG_LEVEL=0)

IF(STATIC_LIBGCC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc")
endif()

if(BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBRARIES OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

function(set_cxx_standard std flag)
  if(${CMAKE_VERSION} VERSION_LESS "3.1.0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
  else()
    set(CMAKE_CXX_STANDARD ${std} PARENT_SCOPE)
  endif()
endfunction(set_cxx_standard)

check_type_size(ptrdiff_t SIZEOF_PTRDIFF_T)
check_type_size(int SIZEOF_INT)

check_cxx_compiler_flag("-std=gnu++14" HAS_GNUPP14)
check_cxx_compiler_flag("-std=gnu++17" HAS_GNUPP17)
check_cxx_compiler_flag("-std=gnu++20" HAS_GNUPP20)
check_cxx_compiler_flag("-std=gnu++23" HAS_GNUPP23)
if(HAS_GNUPP23 OR ${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  set_cxx_standard(23 "-std=gnu++23")
elseif(HAS_GNUPP20)
  set_cxx_standard(20 "-std=gnu++20")
elseif(HAS_GNUPP17)
  set_cxx_standard(17 "-std=gnu++17")
elseif (HAS_GNUPP14)
  set_cxx_standard(14 "-std=gnu++14")
else()
  set_cxx_standard(11 "-std=gnu++11")
endif()

check_symbol_exists(sysinfo "sys/sysinfo.h" HAVE_SYSINFO)

if(CMAKE_BUILD_MARCH)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${CMAKE_BUILD_MARCH}")
endif()

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_HAS_STD_BYTE=0)
  add_definitions(-D_SILENCE_ALL_CXX23_DEPRECATION_WARNINGS)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-implicit-fallthrough -Wreturn-type -Wno-unused -Wno-unused-parameter -Wno-unused-variable -Wno-uninitialized -Wno-deprecated-copy -Wno-unknown-warning-option ")#-g -fsanitize=address -fno-omit-frame-pointer ")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-pragma-clang-attribute -Wno-overloaded-virtual -Wno-missing-braces") #-g -fsanitize=thread -fno-omit-frame-pointer" )
endif()

if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fabi-version=7")
  message("Setting -fabi-version for GCC 4.x")
endif()

set(DISPATCH_OBJECTS
        "src/dp/swipe/banded_3frame_swipe.cpp"
        "src/search/hamming/stage1_2.cpp"
        "src/tools/benchmark.cpp"
        "src/dp/swipe/swipe_wrapper.cpp"
        "src/masking/tantan.cpp"
        "src/dp/scan_diags.cpp"
        "src/dp/ungapped_simd.cpp"
        "src/dp/swipe/anchored_wrapper.cpp"
        "src/dp/score_profile.cpp"
        )

if(EXTRA)
  LIST(APPEND DISPATCH_OBJECTS "src/tools/benchmark_swipe.cpp")
endif()

add_library(arch_generic OBJECT ${DISPATCH_OBJECTS})
target_compile_options(arch_generic PUBLIC -DDISPATCH_ARCH=ARCH_GENERIC -DARCH_ID=0)
if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    target_compile_options(arch_generic PUBLIC /Zc:__cplusplus)
endif()

if(X86)
  add_library(arch_sse4_1 OBJECT ${DISPATCH_OBJECTS})
  add_library(arch_avx2 OBJECT ${DISPATCH_OBJECTS})
  if(WITH_AVX512)
    add_library(arch_avx512 OBJECT ${DISPATCH_OBJECTS})
  endif()
  add_definitions(-DWITH_AVX2)
  add_definitions(-DWITH_SSE4_1)
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    target_compile_options(arch_sse4_1 PUBLIC -DDISPATCH_ARCH=ARCH_SSE4_1 -DARCH_ID=1 -D__SSSE3__ -D__SSE4_1__ -D__POPCNT__ /Zc:__cplusplus)
    target_compile_options(arch_avx2 PUBLIC -DDISPATCH_ARCH=ARCH_AVX2 -DARCH_ID=2 /arch:AVX2 -D__SSSE3__ -D__SSE4_1__ -D__POPCNT__ /Zc:__cplusplus)
    if(WITH_AVX512)
      target_compile_options(arch_avx512 PUBLIC -DDISPATCH_ARCH=ARCH_AVX512 -DARCH_ID=3 /arch:AVX512 -D__SSSE3__ -D__SSE4_1__ -D__POPCNT__)
    endif()
  else()
    target_compile_options(arch_sse4_1 PUBLIC -DDISPATCH_ARCH=ARCH_SSE4_1 -DARCH_ID=1 -mssse3 -mpopcnt -msse4.1)
    target_compile_options(arch_avx2 PUBLIC -DDISPATCH_ARCH=ARCH_AVX2 -DARCH_ID=2 -mssse3 -mpopcnt -msse4.1 -msse4.2 -mavx -mavx2)
    if(WITH_AVX512)
      target_compile_options(arch_avx512 PUBLIC -DDISPATCH_ARCH=ARCH_AVX512 -DARCH_ID=3 -mssse3 -mpopcnt -msse4.1 -msse4.2 -mavx -mavx2 -mavx512f -mavx512bw)
    endif()
  endif()
endif(X86)

# NEON is mandatory on Aarch64
if(AARCH64)
  add_definitions(-DWITH_NEON)
  add_library(arch_neon OBJECT ${DISPATCH_OBJECTS})
  target_compile_options(arch_neon PUBLIC -DDISPATCH_ARCH=ARCH_NEON -DARCH_ID=4 -D__ARM_NEON -D__aarch64__)
endif(AARCH64)

# NEON is optional on Armv7, so we need to check for compiler support,
# and for the <sys/auxv.h> header used for runtime detection.
if(ARM)
check_symbol_exists(getauxval "sys/auxv.h" HAVE_GETAUXVAL)
check_cxx_compiler_flag("-mfpu=neon" HAVE_MFPU_NEON)
if(HAVE_MFPU_NEON)
  message(STATUS "Compiler supports NEON architecture")
    add_definitions(-DWITH_NEON)
    add_definitions(-DHAVE_MFPU_NEON)
    add_library(arch_neon OBJECT ${DISPATCH_OBJECTS})
    target_compile_options(arch_neon PUBLIC -DDISPATCH_ARCH=ARCH_NEON -DARCH_ID=4 -D__ARM_NEON -mfpu=neon)
  endif()
  if(HAVE_GETAUXVAL)
    add_definitions(-DHAVE_GETAUXVAL)
  endif()
endif(ARM)

set(OBJECTS
        src/run/main.cpp
        src/basic/config.cpp
        src/stats/score_matrix.cpp
        src/data/queries.cpp
        src/data/seed_histogram.cpp
        src/output/daa/daa_record.cpp
        src/util/command_line_parser.cpp
        src/util/util.cpp
        src/basic/basic.cpp
        src/basic/hssp.cpp
        src/dp/ungapped_align.cpp
        src/run/tools.cpp
        src/chaining/greedy_align.cpp
        src/output/output_format.cpp
        src/output/join_blocks.cpp
        src/data/frequent_seeds.cpp
        src/align/legacy/query_mapper.cpp
        src/output/blast_tab_format.cpp
        src/output/blast_pairwise_format.cpp
        src/run/double_indexed.cpp
        src/output/sam_format.cpp
        src/align/align.cpp
        src/search/setup.cpp
        src/data/taxonomy.cpp
        src/masking/masking.cpp
        src/data/seed_set.cpp
        src/util/simd.cpp
        src/output/taxon_format.cpp
        src/output/daa/view.cpp
        src/output/output_sink.cpp
        src/output/target_culling.cpp
        src/align/legacy/banded_swipe_pipeline.cpp
        src/util/io/compressed_stream.cpp
        src/util/io/compressed_buffer.cpp
        src/util/io/deserializer.cpp
        src/util/io/file_sink.cpp
        src/util/io/file_source.cpp
        src/util/io/input_file.cpp
        src/util/io/input_stream_buffer.cpp
        src/util/io/output_file.cpp
        src/util/io/output_stream_buffer.cpp
        src/util/io/serializer.cpp
        src/util/io/temp_file.cpp
        src/util/io/text_input_file.cpp
        src/util/io/file.cpp
        src/data/taxon_list.cpp
        src/data/taxonomy_nodes.cpp
        src/lib/murmurhash/MurmurHash3.cpp
        src/search/stage0.cpp
        src/data/seed_array.cpp
        src/output/paf_format.cpp
        src/util/system/system.cpp
        src/util/algo/greedy_vertex_cover.cpp
        src/util/sequence/sequence.cpp
        src/tools/tools.cpp
        src/util/system/getRSS.cpp
        src/lib/tantan/LambdaCalculator.cc
        src/util/string/string.cpp
        src/align/extend.cpp
        src/test/test.cpp
        src/align/ungapped.cpp
        src/align/gapped_score.cpp
        src/align/gapped_final.cpp
        src/align/full_db.cpp
        src/align/culling.cpp
        src/cluster/cluster_registry.cpp
        src/cluster/cascaded/cascaded.cpp
        src/align/output.cpp
        #src/tools/roc.cpp
        src/test/data.cpp
        src/test/test_cases.cpp
        src/chaining/smith_waterman.cpp
        src/output/xml_format.cpp
        src/align/gapped_filter.cpp
        src/util/parallel/filestack.cpp
        src/util/parallel/parallelizer.cpp
        src/util/parallel/multiprocessing.cpp
        src/lib/alp/njn_dynprogprob.cpp
        src/lib/alp/njn_dynprogproblim.cpp
        src/lib/alp/njn_dynprogprobproto.cpp
        src/lib/alp/njn_ioutil.cpp
        src/lib/alp/njn_localmaxstat.cpp
        src/lib/alp/njn_localmaxstatmatrix.cpp
        src/lib/alp/njn_localmaxstatutil.cpp
        src/lib/alp/njn_random.cpp
        src/lib/alp/sls_alignment_evaluer.cpp
        src/lib/alp/sls_alp.cpp
        src/lib/alp/sls_alp_data.cpp
        src/lib/alp/sls_alp_regression.cpp
        src/lib/alp/sls_alp_sim.cpp
        src/lib/alp/sls_basic.cpp
        src/lib/alp/sls_pvalues.cpp
        src/align/global_ranking/global_ranking.cpp
        src/align/global_ranking/extend.cpp
        #src/tools/rocid.cpp
        src/lib/blast/blast_seg.cpp
        src/lib/blast/blast_filter.cpp
        src/lib/blast/nlm_linear_algebra.cpp
        src/stats/stats.cpp
        src/stats/cbs.cpp
        src/stats/comp_based_stats.cpp
        src/stats/hauser_correction.cpp
        src/stats/matrix_adjust.cpp
        src/data/index.cpp
        src/data/dmnd/dmnd.cpp
        src/data/sequence_file.cpp
        src/tools/find_shapes.cpp
        src/data/block/block.cpp
        src/data/block/block_wrapper.cpp
        src/run/config.cpp
        src/data/sequence_set.cpp
        src/align/global_ranking/table.cpp
        src/output/daa/daa_write.cpp
        src/search/seed_complexity.cpp
        src/util/tsv/tsv.cpp
        src/basic/value.cpp
        src/masking/motifs.cpp
        src/align/alt_hsp.cpp
        src/data/fasta/fasta_file.cpp
        src/cluster/output.cpp
        src/cluster/realign.cpp
        src/cluster/reassign.cpp
        src/util/tsv/read_tsv.cpp
        src/tools/greedy_vertex_cover.cpp
        src/cluster/cascaded/recluster.cpp
        src/cluster/helpers.cpp
        src/search/kmer_ranking.cpp
        src/chaining/hamming_ext.cpp
        src/lib/blast/blast_message.cpp
        src/lib/blast/blast_stat.cpp
        src/lib/blast/blastn_score.cpp
        src/lib/blast/ncbi_std.cpp
        src/util/tsv/table.cpp
        src/util/tsv/file.cpp
        src/util/tsv/record.cpp
        src/cluster/cascaded/helpers.cpp
        src/cluster/cascaded/wrapper.cpp
        src/output/daa/merge.cpp
        src/chaining/backtrace.cpp
        src/util/tsv/merge.cpp
        src/util/tsv/join.cpp
        src/dp/scalar/smith_waterman.cpp
        src/tools/tsv.cpp
        src/cluster/external/external.cpp
        src/cluster/external/align.cpp
        src/cluster/external/cluster.cpp
        src/cluster/external/output.cpp
        #src/util/algo/mcl.cpp
        #src/stats/blast/matrix_adjust.cpp
        src/stats/blast/ncbi.cpp
        src/search/hit_buffer.cpp
        src/data/blastdb/phr.cpp
        src/data/blastdb/psq.cpp
        src/data/blastdb/pin.cpp
        src/data/blastdb/pal.cpp
        src/data/blastdb/asn1.cpp
        src/data/blastdb/blastdb.cpp
)

if(WITH_DNA)
  list(APPEND OBJECTS
          src/contrib/dna/smith_watermann.cpp
          src/contrib/dna/dna_index.cpp
          src/contrib/dna/seed_set_dna.cpp
          src/contrib/dna/extension.cpp
          src/contrib/dna/chain.cpp
          src/contrib/dna/extension_chain.cpp
          src/contrib/dna/alignment.cpp
          src/contrib/dna/extension_seed_matches.cpp
          src/contrib/dna/smith_watermann.cpp
          src/lib/ksw2/ksw2_extz2_sse.c
          src/lib/ksw2/ksw2_extz.c
          src/lib/WFA2-lib.diamond/bindings/cpp/WFAligner.cpp
          src/contrib/dna/build_score.cpp
  )
endif()

if(WITH_MCL)
  list(APPEND OBJECTS src/contrib/mcl/mcl.cpp
                      src/contrib/mcl/clustering_variables.cpp
                      src/contrib/mcl/clustering_format.cpp)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

if(WITH_ZSTD)
  set(ZSTD_OBJ "src/util/io/zstd_stream.cpp")
endif()

if(X86)
  if(WITH_AVX512)
    add_executable(diamond $<TARGET_OBJECTS:arch_generic> $<TARGET_OBJECTS:arch_sse4_1> $<TARGET_OBJECTS:arch_avx2> $<TARGET_OBJECTS:arch_avx512> ${OBJECTS} ${ZSTD_OBJ})
  else()
    add_executable(diamond $<TARGET_OBJECTS:arch_generic> $<TARGET_OBJECTS:arch_sse4_1> $<TARGET_OBJECTS:arch_avx2> ${OBJECTS} ${ZSTD_OBJ})
  endif()
elseif(ARM OR AARCH64)
  if(HAVE_MFPU_NEON)
     message(STATUS "Using NEON architecture")
     add_executable(diamond $<TARGET_OBJECTS:arch_generic> $<TARGET_OBJECTS:arch_neon> ${OBJECTS} ${ZSTD_OBJ})
  else()
     add_executable(diamond $<TARGET_OBJECTS:arch_generic> ${OBJECTS} ${ZSTD_OBJ})
  endif()
else()
  add_executable(diamond $<TARGET_OBJECTS:arch_generic> ${OBJECTS} ${ZSTD_OBJ})
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  target_compile_options(diamond PUBLIC /Zc:__cplusplus)
endif()

include_directories("${CMAKE_SOURCE_DIR}/src/" "${ZLIB_INCLUDE_DIR}" "${CMAKE_SOURCE_DIR}/src/lib")
#find_package(SQLite3)
#target_link_libraries(diamond SQLite::SQLite3 ${CMAKE_DL_LIBS})
find_path(${SQLite3_INCLUDE_DIR} NAMES sqlite3.h PATHS ${SQLite3_INCLUDE_DIR})
find_library(SQLite3_LIBRARY NAMES libsqlite3.a libsqlite3.so libsqlite3.so.0 sqlite3 HINTS /opt/homebrew/opt/sqlite/lib /usr/local/opt/sqlite/lib)
message(STATUS "SQLite3 Library: ${SQLite3_LIBRARY}")
target_include_directories(diamond PRIVATE "${SQLite3_INCLUDE_DIR}")
target_link_libraries(diamond ${SQLite3_LIBRARY} ${CMAKE_DL_LIBS})

if(WITH_ZSTD)
  find_path(${ZSTD_INCLUDE_DIR} NAMES zstd.h PATHS ${ZSTD_INCLUDE_DIR})
  find_library(ZSTD_LIBRARY NAMES libzstd.a libzstd_static)
  target_include_directories(diamond PRIVATE "${ZSTD_INCLUDE_DIR}")
  target_link_libraries(diamond ${ZSTD_LIBRARY})
  add_definitions(-DWITH_ZSTD)
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  find_library(BCRYPT_LIBRARY BCrypt)
  find_library(WS2_32_LIBRARY Ws2_32)
  target_link_libraries(diamond ${BCRYPT_LIBRARY} ${WS2_32_LIBRARY})
endif()

target_link_libraries(diamond ${ZLIB_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS diamond DESTINATION bin)

enable_testing()
SET(TD ${CMAKE_SOURCE_DIR}/src/test)
SET(SP -DTEST_DIR=${CMAKE_SOURCE_DIR}/src/test -P ${CMAKE_SOURCE_DIR}/src/test/test.cmake)
add_test(NAME blastp COMMAND ${CMAKE_COMMAND} -DNAME=blastp "-DARGS=blastp -q ${TD}/1.faa -d ${TD}/2.faa -p1" ${SP})
add_test(NAME blastp-mid-sens COMMAND ${CMAKE_COMMAND} -DNAME=blastp-mid-sens "-DARGS=blastp -q ${TD}/3.faa -d ${TD}/4.faa --mid-sensitive -p1" ${SP})
add_test(NAME blastp-f0 COMMAND ${CMAKE_COMMAND} -DNAME=blastp-f0 "-DARGS=blastp -q ${TD}/1.faa -d ${TD}/2.faa -f0 -p1" ${SP})
add_test(NAME diamond COMMAND diamond test)
