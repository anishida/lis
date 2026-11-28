noinst_LTLIBRARIES = libsolver.la
libsolver_la_SOURCES = lis_solver.c lis_solver_nearnullspace.c \
lis_solver_nullspace_candidates.c lis_solver_component_nullspace.c \
lis_solver_nullspace_lift.c lis_solver_nullspace_query.c \
lis_solver_capabilities.c lis_solver_options.c \
lis_solver_bicg.c \
lis_solver_bicgsafe.c \
lis_solver_bicgstab.c \
lis_solver_bicgstabl.c \
lis_solver_cg.c \
lis_solver_cgs.c \
lis_solver_gmres.c \
lis_solver_gcrodr.c \
lis_solver_fgcrodr.c \
lis_solver_gpbicg.c \
lis_solver_gs.c \
lis_solver_idrs.c \
lis_solver_jacobi.c \
lis_solver_minres.c \
lis_solver_orthomin.c \
lis_solver_qmr.c \
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
LTCOMPILE = $(LIBTOOL) --mode=compile --tag=CC $(CC) $(AMDEFS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)

# Native RCI generation is local to this build directory.
if ENABLE_NATIVE_RCI
RCI_NATIVE_INPUTS = $(top_srcdir)/src/solver/lis_solver.c \
  $(top_srcdir)/src/solver/lis_solver_cg.c \
  $(top_srcdir)/src/solver/lis_solver_bicgstab.c \
  $(top_srcdir)/src/solver/lis_solver_gmres.c \
  $(top_srcdir)/src/solver/lis_solver_gcrodr.c \
  $(top_srcdir)/src/solver/lis_solver_fgcrodr.c \
  $(top_srcdir)/src/esolver/lis_esolver.c \
  $(top_srcdir)/src/esolver/lis_esolver_block.c \
  $(top_srcdir)/src/esolver/lis_esolver_feast.c \
  $(top_srcdir)/src/esolver/lis_krylov.c \
  $(top_srcdir)/qualification/rci-restore/generate_krylov_continuation.js \
  $(top_srcdir)/qualification/rci-restore/generate_linear_continuation.py \
  $(top_srcdir)/qualification/rci-restore/krylov_continuation.inc \
  $(top_srcdir)/qualification/rci-restore/lis_rci_augmented_restore.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_augmented_restore.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_feast_restore.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_feast_restore.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_krylov_restore.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_krylov_restore.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_feast.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_feast.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_feast_lifecycle.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_feast_lifecycle.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_feast_projector.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_feast_projector.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_fgcro_restore.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_fgcro_restore.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_gcro_restore.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_gcro_restore.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_lifecycle.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_lifecycle.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_linear_restore.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_linear_restore.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_public.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_restore.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_restore.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_session.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_native_session.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_recycle_restore.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_recycle_restore.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_restore.c \
  $(top_srcdir)/qualification/rci-restore/lis_rci_restore.h \
  $(top_srcdir)/qualification/rci-restore/lis_rci_split_restore.c \
  $(top_srcdir)/qualification/rci-restore/prepare_native_installed_api.py \
  $(top_srcdir)/qualification/rci-restore/qualify_all_native.py \
  $(top_srcdir)/qualification/rci-restore/qualify_feast_native.py \
  $(top_srcdir)/qualification/rci-restore/qualify_feast_projector_native.py \
  $(top_srcdir)/qualification/rci-restore/qualify_fgcro_native.py \
  $(top_srcdir)/qualification/rci-restore/qualify_gcro_native.py \
  $(top_srcdir)/qualification/rci-restore/qualify_krylov_native.py \
  $(top_srcdir)/qualification/rci-restore/qualify_linear_native.py \
  $(top_srcdir)/qualification/rci-restore/qualify_native_feast_lifecycle.py \
  $(top_srcdir)/qualification/rci-restore/qualify_native_lifecycle.py \
  $(top_srcdir)/qualification/rci-restore/build_native_tree.py \
  $(top_srcdir)/include/lis_rci_native.h
rci-native.stamp: $(RCI_NATIVE_INPUTS) Makefile
	python3 "$(top_srcdir)/qualification/rci-restore/build_native_tree.py" "$(top_srcdir)" "$(abs_builddir)/rci-generated" "$(top_builddir)/include/lis_config.h"
	touch $@
lis_solver.lo: rci-native.stamp $(top_builddir)/include/lis_config.h
	$(LTCOMPILE) -I$(top_builddir)/include -c -o $@ rci-generated/src/solver/lis_solver.c
lis_solver_cg.lo: rci-native.stamp $(top_builddir)/include/lis_config.h
	$(LTCOMPILE) -I$(top_builddir)/include -c -o $@ rci-generated/src/solver/lis_solver_cg.c
lis_solver_bicgstab.lo: rci-native.stamp $(top_builddir)/include/lis_config.h
	$(LTCOMPILE) -I$(top_builddir)/include -c -o $@ rci-generated/src/solver/lis_solver_bicgstab.c
lis_solver_gmres.lo: rci-native.stamp $(top_builddir)/include/lis_config.h
	$(LTCOMPILE) -I$(top_builddir)/include -c -o $@ rci-generated/src/solver/lis_solver_gmres.c
lis_solver_gcrodr.lo: rci-native.stamp $(top_builddir)/include/lis_config.h
	$(LTCOMPILE) -I$(top_builddir)/include -c -o $@ rci-generated/src/solver/lis_solver_gcrodr.c
lis_solver_fgcrodr.lo: rci-native.stamp $(top_builddir)/include/lis_config.h
	$(LTCOMPILE) -I$(top_builddir)/include -c -o $@ rci-generated/src/solver/lis_solver_fgcrodr.c
endif

clean-local:
	rm -rf rci-generated rci-native.stamp

# Recompile after changing configure options, including disabling RCI.
lis_solver.lo lis_solver_cg.lo lis_solver_bicgstab.lo lis_solver_gmres.lo lis_solver_gcrodr.lo lis_solver_fgcrodr.lo: Makefile $(top_builddir)/include/lis_config.h
