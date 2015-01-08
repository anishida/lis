noinst_LTLIBRARIES = libsystem.la
libsystem_la_SOURCES = lis_error.c \
lis_hash.c \
lis_init.c \
lis_input.c \
lis_input_hb.c \
lis_input_mm.c \
lis_memory.c \
lis_output.c \
lis_output_mm.c \
lis_sort.c \
lis_time.c \
mt19937ar.c


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
