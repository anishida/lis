noinst_LTLIBRARIES = libmatvec.la
libmatvec_la_SOURCES = lis_matvec.c     \
lis_matvec_bsc.c \
lis_matvec_bsr.c \
lis_matvec_csc.c \
lis_matvec_coo.c \
lis_matvec_csr.c \
lis_matvec_dia.c \
lis_matvec_dns.c \
lis_matvec_ell.c \
lis_matvec_jad.c \
lis_matvec_msr.c \
lis_matvec_vbr.c


AM_CFLAGS = -I$(top_srcdir)/include
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
