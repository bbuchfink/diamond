cmake_minimum_required (VERSION 2.6)
project (DIAMOND)

option(BUILD_STATIC "BUILD_STATIC" OFF)

if(BUILD_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBRARIES OFF)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

find_package(Boost
  1.53.0 REQUIRED
  COMPONENTS program_options timer iostreams thread
  system timer chrono)

find_package(ZLIB REQUIRED)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-uninitialized")

add_subdirectory("src/algo/blast/core")
include_directories(
  "${CMAKE_SOURCE_DIR}/src"
  "${Boost_INCLUDE_DIR}"
  "${ZLIB_INCLUDE_DIR}")

add_executable(diamond src/main.cpp
  src/basic/options.cpp
  src/util/tinythread.cpp)

target_link_libraries(diamond blast_core ${Boost_LIBRARIES} ${ZLIB_LIBRARY})

install(TARGETS diamond DESTINATION bin)
