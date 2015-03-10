ACLOCAL_AMFLAGS = -I m4
#EXTRA_DIST = autogen.sh
bin_PROGRAMS = diamond
diamond_SOURCES = algo/blast/core/blast_encoding.c algo/blast/core/blast_stat.c algo/blast/core/blast_filter.c algo/blast/core/blast_util.c algo/blast/core/blast_message.c algo/blast/core/ncbi_erf.c algo/blast/core/blast_options.c algo/blast/core/ncbi_math.c algo/blast/core/blast_program.c algo/blast/core/ncbi_std.c algo/blast/core/blast_psi_priv.c algo/blast/core/raw_scoremat.c algo/blast/core/blast_query_info.c algo/blast/core/blast_seg.c main.cpp basic/options.cpp
diamond_CPPFLAGS = -DNDEBUG $(BOOST_CPPFLAGS)
diamond_LDFLAGS = $(BOOST_THREAD_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_TIMER_LDFLAGS) $(BOOST_CHRONO_LDFLAGS) $(BOOST_SYSTEM_LDFLAGS) $(BOOST_IOSTREAMS_LDFLAGS) -all-static -fopenmp
diamond_LDADD = $(BOOST_THREAD_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_TIMER_LIBS) $(BOOST_CHRONO_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) -lrt -lz
