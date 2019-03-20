AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ -I$(srcdir)/../lib
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) $(AM_LDFLAGS) ../lib/libbiomcmc.la

LIST_OF_TEST_PROGS= check_unit
# duplicate, since we want all to be compiled only with 'make check'
TESTS = $(LIST_OF_TEST_PROGS)           # list of test programs 
check_PROGRAMS = $(LIST_OF_TEST_PROGS)  # list of programs to be compiled only with 'make check'

check_unit_SOURCES = check_unit.c
