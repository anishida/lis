noinst_LTLIBRARIES = libesolver.la
libesolver_la_SOURCES = lis_esolver.c lis_gesolve_matrix_session.c \
lis_esolver_config.c \
lis_esolver_selection.c lis_esolver_selection_copy.c \
lis_esolver_block.c \
lis_esolver_ai.c \
lis_esolver_cg.c \
lis_esolver_ii.c \
lis_esolver_li.c \
lis_esolver_pi.c \
lis_esolver_rqi.c \
lis_esolver_feast.c \
lis_esolver_si.c	

libesolver_la_SOURCES += lis_krylov.c lis_esolver_krylov.c

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
lis_esolver.lo: rci-native.stamp $(top_builddir)/include/lis_config.h
	$(LTCOMPILE) -I$(top_builddir)/include -c -o $@ rci-generated/src/esolver/lis_esolver.c
lis_esolver_block.lo: rci-native.stamp $(top_builddir)/include/lis_config.h
	$(LTCOMPILE) -I$(top_builddir)/include -c -o $@ rci-generated/src/esolver/lis_esolver_block.c
lis_esolver_feast.lo: rci-native.stamp $(top_builddir)/include/lis_config.h
	$(LTCOMPILE) -I$(top_builddir)/include -c -o $@ rci-generated/src/esolver/lis_esolver_feast.c
lis_krylov.lo: rci-native.stamp $(top_builddir)/include/lis_config.h
	$(LTCOMPILE) -I$(top_builddir)/include -c -o $@ rci-generated/src/esolver/lis_krylov.c
endif

clean-local:
	rm -rf rci-generated rci-native.stamp

# Recompile after changing configure options, including disabling RCI.
lis_esolver.lo lis_esolver_block.lo lis_esolver_feast.lo lis_krylov.lo: Makefile $(top_builddir)/include/lis_config.h
