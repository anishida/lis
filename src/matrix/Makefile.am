noinst_LTLIBRARIES = libmatrix.la
libmatrix_la_SOURCES = \
lis_matrix.c               \
lis_matrix_bsc.c           \
lis_matrix_bsr.c           \
lis_matrix_csc.c           \
lis_matrix_coo.c           \
lis_matrix_csr.c           \
lis_matrix_dia.c           \
lis_matrix_diag.c          \
lis_matrix_dns.c           \
lis_matrix_ell.c           \
lis_matrix_ilu.c           \
lis_matrix_jad.c           \
lis_matrix_mpi.c           \
lis_matrix_msr.c           \
lis_matrix_ops.c           \
lis_matrix_rco.c           \
lis_matrix_vbr.c           

AM_CFLAGS = -I$(top_srcdir)/include $(ILIBS)
CLEANFILES = *.il

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
