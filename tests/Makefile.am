AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ -I$(srcdir)/../lib -DTEST_FILE_DIR='"$(srcdir)/files/"' 
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@  @CHECK_CFLAGS@ @ZLIB_CFLAGS@  @LZMA_CFLAGS@
LDADD =  ../lib/libbiomcmc_static.la $(GTKDEPS_LIBS) $(AM_LDFLAGS) @CHECK_LIBS@ @ZLIB_LIBS@  @LZMA_LIBS@

EXTRA_DIST = files # directory with fasta etc files (accessed with #define TEST_FILE_DIR above)
# we use the list twice below, since we want all to be compiled only with 'make check'
LIST_OF_TEST_PROGS= check_unit check_topology debug_topology debug_rng debug_gff3 debug_compression debug_hash 

TESTS = $(LIST_OF_TEST_PROGS)           # list of test programs 
check_PROGRAMS = $(LIST_OF_TEST_PROGS)  # list of programs to be compiled only with 'make check' (like noinst_PROGRAMS)

#check_minhash_SOURCES = check_minhash.c
#check_suffix_tree_SOURCES = check_suffix_tree.c
check_unit_SOURCES = check_unit.c # ../lib/config.h   ## config.h must be mentioned at least once 
check_topology_SOURCES = check_topology.c
# not using libcheck, not actual tests
debug_topology_SOURCES = debug_topology.c
debug_rng_SOURCES = debug_rng.c
debug_gff3_SOURCES = debug_gff3.c
debug_compression_SOURCES = debug_compression.c
debug_hash_SOURCES = debug_hash.c
