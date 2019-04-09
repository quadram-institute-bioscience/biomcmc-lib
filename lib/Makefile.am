AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) $(AM_LDFLAGS)

# ideally I would hide reconciliation.h etc used only inside genetree.c 
common_headers = hashtable.h lowlevel.h random_number_gen.h random_number.h nexus_common.h \
                 bipartition.h prob_distribution.h empirical_frequency.h biomcmc.h argtable3.h \
                 distance_matrix.h alignment.h topology_common.h parsimony.h genetree.h \
								 reconciliation.h splitset_distances.h 

common_src     = hashtable.c lowlevel.c random_number_gen.c random_number.c nexus_common.c \
                 bipartition.c prob_distribution.c empirical_frequency.c argtable3.c \
                 distance_matrix.c alignment.c topology_common.c parsimony.c \
								 reconciliation.c splitset_distances.c genetree.c 


noinst_LTLIBRARIES = libbiomcmc.la # convenience lib (not installed globally)
#lib_LTLIBRARIES = libbiomcmc.la  # libtools 
#otherincludedir = $(includedir)/biomcmc
#otherinclude_HEADERS = config.h $(common_headers) # if headers are here (=global) should not be on SOURCES (=local)
#libbiomcmc_la_SOURCES = $(common_src) $(common_headers) config.h
libbiomcmc_la_SOURCES = $(common_src) config.h
