cmake_minimum_required (VERSION 2.6)
project (DIAMOND)

find_package(Boost
  1.53.0 REQUIRED
  COMPONENTS program_options timer iostreams thread
  system timer chrono)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-uninitialized")

add_subdirectory("src/algo/blast/core")
include_directories(
  "${CMAKE_SOURCE_DIR}/src"
  "${Boost_INCLUDE_DIR}")

add_executable(diamond src/main.cpp
  src/basic/options.cpp
  src/util/tinythread.cpp)

target_link_libraries(diamond blast_core ${Boost_LIBRARIES})

install(TARGETS diamond DESTINATION bin)
