#!/bin/sh

# Sample shell script to draw residual histories using gnuplot
# Make executable files and run the following command:
#
# > ./gunplot1b.sh
#


# Check that the utility is present.

type gnuplot >/dev/null 2>&1 || { echo >&2 "You need gnuplot to run the script. Aborting."; exit 1; }


# Run test programs.

srcdir=../test
$srcdir/test1 $srcdir/testmat.mtx 1 /dev/null ILU0-BiCG -i bicg -p ilu
$srcdir/test1 $srcdir/testmat.mtx 1 /dev/null ILU0-CGS -i cgs -p ilu
$srcdir/test1 $srcdir/testmat.mtx 1 /dev/null ILU0-BiCGSTAB -i bicgstab -p ilu


# Draw residual histories using gnuplot.

gnuplot -e "filenames='ILU0-BiCG ILU0-CGS ILU0-BiCGSTAB'" gnuplot1b.plt


rm ILU0-*

