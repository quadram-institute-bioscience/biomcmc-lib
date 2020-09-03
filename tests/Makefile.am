AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ -I$(srcdir)/../lib -DTEST_FILE_DIR='"$(srcdir)/files/"' 
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@  @CHECK_CFLAGS@ @ZLIB_CFLAGS@ 
LDADD = $(GTKDEPS_LIBS) $(AM_LDFLAGS) @CHECK_LIBS@ @ZLIB_LIBS@ ../lib/libbiomcmc.la 

EXTRA_DIST = files # directory with fasta etc files (accessed with #define TEST_FILE_DIR above)
# we use the list twice below, since we want all to be compiled only with 'make check'
LIST_OF_TEST_PROGS= check_unit check_topology debug_topology debug_rng 

TESTS = $(LIST_OF_TEST_PROGS)           # list of test programs 
check_PROGRAMS = $(LIST_OF_TEST_PROGS)  # list of programs to be compiled only with 'make check' (like noinst_PROGRAMS)

#check_minhash_SOURCES = check_minhash.c
#check_suffix_tree_SOURCES = check_suffix_tree.c
check_unit_SOURCES = check_unit.c
check_topology_SOURCES = check_topology.c
# not using libcheck, not actual tests
debug_topology_SOURCES = debug_topology.c
debug_rng_SOURCES = debug_rng.c
