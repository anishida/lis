#!/bin/sh

# Sample shell script to draw eigenvectors using gnuplot
# Make executable files and run the following command:
#
# > ./gunplot5.sh
#


# Check that the utility is present.

type gnuplot >/dev/null 2>&1 || { echo >&2 "You need gnuplot to run the script. Aborting."; exit 1; }


# Run test programs.

subspace=3
srcdir=../test
$srcdir/etest5 $srcdir/testmat.mtx /dev/null evectors.mtx /dev/null /dev/null -e si -ie ii -ss $subspace


# Draw eigenvectors.

filename=evectors.mtx
while read M N L B X; do
    echo $M | grep -v '^%.*' > /dev/null
    if [ $? -eq 0 ]; then
	size=$M
	break
    fi
done < $filename

gnuplot -e "filename='$filename'; size='$size'; subspace='$subspace'" gnuplot5.plt

rm evectors.mtx
