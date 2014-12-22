noinst_LTLIBRARIES = libsolver.la
libsolver_la_SOURCES = lis_solver.c           \
lis_solver_bicg.c      \
lis_solver_bicgsafe.c  \
lis_solver_bicgstab.c  \
lis_solver_bicgstabl.c \
lis_solver_cg.c        \
lis_solver_cgs.c       \
lis_solver_gmres.c     \
lis_solver_gpbicg.c    \
lis_solver_gs.c        \
lis_solver_idrs.c      \
lis_solver_jacobi.c    \
lis_solver_minres.c    \
lis_solver_orthomin.c  \
lis_solver_qmr.c       \
lis_solver_sor.c       


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
