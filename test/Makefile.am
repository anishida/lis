LDADD = -L$(top_builddir)/src -llis@LIBSUFFIX@
AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/include
AM_FFLAGS = -I$(top_builddir) -I$(top_srcdir)/include
AM_LDFLAGS = $(AMDEFS)

EXTRA_DIST = defs test.sh testmat.mtx testmat2.mtx testmat3.mtx testmat4.mtx testmat0.mtx testvec0.mtx

testdir = ${prefix}/share/examples/lis

test_DATA = testmat.mtx testmat2.mtx testmat3.mtx testmat4.mtx testmat0.mtx testvec0.mtx

test_SCRIPTS = defs test.sh

test_PROGRAMS = etest1 etest2 etest3 etest4 etest5 etest5b etest6 etest7 getest1 getest5 getest5b spmvtest1 spmvtest2 spmvtest2b spmvtest3 spmvtest3b spmvtest4 spmvtest5 test1 test2 test2b test3 test3b test3c test4 test5 test6 test7
if ENABLE_FORTRAN
  test_PROGRAMS += etest1f etest4f test1f test4f test7f
endif
if ENABLE_F90
  test_PROGRAMS += test2f test6f 
if ENABLE_SAAMG
if ENABLE_MPI
  test_PROGRAMS += test8f
endif
endif
endif

etest1_SOURCES = etest1.c
etest2_SOURCES = etest2.c
etest3_SOURCES = etest3.c
etest4_SOURCES = etest4.c
etest5_SOURCES = etest5.c
etest5b_SOURCES = etest5b.c
etest6_SOURCES = etest6.c
etest7_SOURCES = etest7.c
getest1_SOURCES = getest1.c
getest5_SOURCES = getest5.c
getest5b_SOURCES = getest5b.c
spmvtest1_SOURCES = spmvtest1.c
spmvtest2_SOURCES = spmvtest2.c
spmvtest2b_SOURCES = spmvtest2b.c
spmvtest3_SOURCES = spmvtest3.c
spmvtest3b_SOURCES = spmvtest3b.c
spmvtest4_SOURCES = spmvtest4.c
spmvtest5_SOURCES = spmvtest5.c
test1_SOURCES = test1.c
test2_SOURCES = test2.c
test2b_SOURCES = test2b.c
test3_SOURCES = test3.c
test3b_SOURCES = test3b.c
test3c_SOURCES = test3c.c
test4_SOURCES = test4.c
test5_SOURCES = test5.c
test6_SOURCES = test6.c
test7_SOURCES = test7.c
if ENABLE_FORTRAN
  etest1f_SOURCES = etest1f.F
  etest4f_SOURCES = etest4f.F
  test1f_SOURCES = test1f.F
  test4f_SOURCES = test4f.F
  test7f_SOURCES = test7f.F
endif
if ENABLE_F90
  test2f_SOURCES = test2f.F90
  test6f_SOURCES = test6f.F90
  test8f_SOURCES = test8f.F90
endif

bin_PROGRAMS = esolve esolver gesolve gesolver hpcg_kernel hpcg_spmvtest lsolve spmvtest1 spmvtest2 spmvtest2b spmvtest3 spmvtest3b spmvtest4 spmvtest5

lsolve_SOURCES = test1.c
esolve_SOURCES = etest5.c
esolver_SOURCES = etest5b.c
gesolve_SOURCES = getest5.c
gesolver_SOURCES = getest5b.c
hpcg_kernel_SOURCES = test3b.c
hpcg_spmvtest_SOURCES = spmvtest3b.c

PPF77COMPILE = $(F77) $(AMDEFS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FFLAGS) $(FFLAGS)
LTPPF77COMPILE = $(LIBTOOL) --mode=compile --tag=F77 $(F77) $(AMDEFS) $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_FFLAGS) $(FFLAGS)
COMPILE = $(CC) $(AMDEFS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
	$(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
LTCOMPILE = $(LIBTOOL) --mode=compile --tag=CC $(CC) $(AMDEFS) $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)

if ENABLE_SAAMG
  LINK = $(LIBTOOL) --mode=link $(FC) \
	$(AM_LDFLAGS) $(LDFLAGS) $(FCLDFLAGS) -o $@ $(LIBS) $(CLIBS)
  F77LINK = $(LIBTOOL) --mode=link $(FC) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@ $(LIBS) $(CLIBS)
else
if ENABLE_FORTRAN
  LINK = $(LIBTOOL) --mode=link $(FC) \
	$(AM_LDFLAGS) $(LDFLAGS) $(FCLDFLAGS) -o $@ $(LIBS) $(CLIBS)
  F77LINK = $(LIBTOOL) --mode=link $(F77) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@ $(LIBS) $(CLIBS)
else
  LINK = $(LIBTOOL) --mode=link $(CC) $(AM_LDFLAGS) $(LDFLAGS) -o $@ $(LIBS)
endif
endif

CLEANFILES = *.il


TESTS_ENVIRONMENT = top_builddir=$(top_builddir) testdir=$(testdir) \
                    enable_f90=$(enable_f90) \
                    enable_fortran=$(enable_fortran) \
                    enable_mpi=$(enable_mpi) \
                    enable_omp=$(enable_omp) \
                    enable_quad=$(enable_quad) \
                    enable_longdouble=$(enable_longdouble) \
                    enable_complex=$(enable_complex) \
                    enable_saamg=$(enable_saamg) \
                    MPIRUN=$(MPIRUN) MPINP=$(MPINP)

TESTS = test.sh

installcheck: 
	@if test ! -f $(testdir)/test.sh; then \
	    echo '=== Test programs not installed'; \
	    exit; \
	fi; \
	$(TESTS_ENVIRONMENT) $(testdir)/test.sh

