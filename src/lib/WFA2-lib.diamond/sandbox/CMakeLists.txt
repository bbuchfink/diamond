cmake_minimum_required(VERSION 3.23)
project(ksw2_Sandbox)

set(CMAKE_CXX_STANDARD 14)

add_executable(ksw2_Sandbox main.cpp
                            ksw2_extz2_sse.c)
