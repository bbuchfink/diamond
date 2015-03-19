dnl Process this file with autoconf to produce a configure script.
AC_INIT(diamond, 0.7.4, buchfink@gmail.com)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR(./main.cpp)

m4_pattern_allow(BOOST_FILESYSTEM_LIBS)
AC_PROG_LIBTOOL
BOOST_REQUIRE([1.53.0])
BOOST_SYSTEM([s])
BOOST_THREAD([s])
BOOST_PROGRAM_OPTIONS([s])
BOOST_CHRONO([s])
BOOST_TIMER([s])
BOOST_IOSTREAMS([s])

AX_CHECK_COMPILE_FLAG([-march=native], [CXXFLAGS="$CXXFLAGS -march=native"])

AM_INIT_AUTOMAKE

AC_PROG_CXX

AC_OUTPUT(Makefile)
