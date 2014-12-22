noinst_LTLIBRARIES = libvector.la
libvector_la_SOURCES = lis_vector.c      \
lis_vector_ops.c  \
lis_vector_opv.c


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
