#!/bin/sh

if test -z "$srcdir"; then
  srcdir=`echo "$0" | sed 's,[^/]*$,,'`
  test "$srcdir" = "$0" && srcdir=.
  test -z "$srcdir" && srcdir=.
  test "${VERBOSE+set}" != set && VERBOSE=1
fi
  VERBOSE=1
. $srcdir/defs

OMP_NUM_THREADS=1
export OMP_NUM_THREADS
if test "$enable_mpi" = "yes"; then
  if test -z "$MPIRUN"; then
    MPIRUN="mpirun -np 2"
  else
    MPIRUN="$MPIRUN $MPINP 2"
  fi
else
  MPIRUN=""
fi
if test "$enable_omp" = "yes"; then
  OMP_NUM_THREADS=2
  export OMP_NUM_THREADS
fi

echo ' '
echo 'checking linear solvers...'
$MPIRUN $srcdir/test1 $srcdir/testmat.mtx 0 /dev/null /dev/null

echo 'checking eigensolvers...'
$MPIRUN $srcdir/etest1 $srcdir/testmat.mtx /dev/null /dev/null

if test "$enable_fortran" = "yes" || test "$enable_f90" = "yes"; then
  echo 'checking Fortran interface...'
  $MPIRUN $srcdir/test4f
  echo ' '
fi

if test "$enable_quad" = "yes"; then
  echo 'checking double-double precision operations...'
  $MPIRUN $srcdir/test5 200 2.0 -f double
  $MPIRUN $srcdir/test5 200 2.0 -f quad
fi

if test "$enable_saamg" = "yes"; then
  echo 'checking SAAMG preconditioner...'
  $MPIRUN $srcdir/test2 10 10 1 /dev/null /dev/null -i cg -p saamg
fi

