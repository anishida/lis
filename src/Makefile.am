SUBDIRS = array esolver matrix matvec precision precon solver system vector 
if ENABLE_SAAMG
  SUBDIRS += fortran
  SUBDIRS += fortran/amg
else
if ENABLE_FORTRAN
  SUBDIRS += fortran
endif
endif

lib_LTLIBRARIES  = liblis@LIBSUFFIX@.la
liblis@LIBSUFFIX@_la_LIBADD = array/libarray.la \
                   esolver/libesolver.la \
                   matrix/libmatrix.la \
                   matvec/libmatvec.la \
                   precision/libprecision.la \
                   precon/libprecon.la \
                   solver/libsolver.la \
                   system/libsystem.la \
                   vector/libvector.la
if ENABLE_SAAMG
  liblis@LIBSUFFIX@_la_LIBADD += fortran/amg/libsaamg.la
  liblis@LIBSUFFIX@_la_LIBADD += fortran/libfortran.la
else
if ENABLE_FORTRAN
  liblis@LIBSUFFIX@_la_LIBADD += fortran/libfortran.la
endif
endif
liblis@LIBSUFFIX@_la_SOURCES = 
