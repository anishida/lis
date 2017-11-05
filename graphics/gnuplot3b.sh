#!/bin/sh

# Sample shell script to draw Ritz value distribution on complex plane
# using gnuplot
# Make executable files enabling complex arithmetic and 
# run the following command:
#
# > ./gunplot3b.sh
#


# Check that the utility is present.

type gnuplot >/dev/null 2>&1 || { echo >&2 "You need gnuplot to run the script. Aborting."; exit 1; }


# Run test programs.

srcdir=../test
$srcdir/etest5b $srcdir/testmat2.mtx rvalues.mtx -e ai -ss 20


# Draw Ritz value distribution.

filename=rvalues.mtx
while read M N L B X; do
    echo $M | grep -v '^%.*' > /dev/null
    if [ $? -eq 0 ]; then
	size=$M
	break
    fi
done < $filename

gnuplot -e "filename='$filename'; size='$size'" gnuplot3b.plt

rm rvalues.mtx
