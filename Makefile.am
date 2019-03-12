AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) $(AM_LDFLAGS)

common_headers = hashtable.h lowlevel.h random_number_gen.h random_number.h nexus_common.h \
                 bipartition.h prob_distribution.h empirical_frequency.h super_sptree.h

common_src     = hashtable.c lowlevel.c random_number_gen.c random_number.c nexus_common.c \
                 bipartition.c prob_distribution.c empirical_frequency.c


lib_LTLIBRARIES = libbiomcmc.la # in automake doc, it's lib_LIBRARIES = libbiomcmc.a
#otherincludedir = $(includedir)/biomcmc
#otherinclude_HEADERS = config.h $(common_headers)
libbiomcmc_la_SOURCES = $(common_src) $(common_headers) config.h
