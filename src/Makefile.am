SUBDIRS = vector matrix matvec array precon solver esolver system precision
if ENABLE_SAAMG
  SUBDIRS += fortran
  SUBDIRS += fortran/amg
else
if ENABLE_FORTRAN
  SUBDIRS += fortran
endif
endif

lib_LTLIBRARIES  = liblis@LIBSUFFIX@.la
liblis@LIBSUFFIX@_la_LIBADD = vector/libvector.la \
                   matrix/libmatrix.la \
                   matvec/libmatvec.la \
                   array/libarray.la \
                   precon/libprecon.la \
                   solver/libsolver.la \
                   esolver/libesolver.la \
                   system/libsystem.la \
                   precision/libprecision.la 
if ENABLE_SAAMG
  liblis@LIBSUFFIX@_la_LIBADD += fortran/libfortran.la
  liblis@LIBSUFFIX@_la_LIBADD += fortran/amg/libsaamg.la
else
if ENABLE_FORTRAN
  liblis@LIBSUFFIX@_la_LIBADD += fortran/libfortran.la
endif
endif
liblis@LIBSUFFIX@_la_SOURCES = 
