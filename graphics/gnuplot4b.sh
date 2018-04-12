#!/bin/sh

# Sample shell script to draw distribution of Ritz values on complex plane
# using gnuplot
# Make executable files enabling complex arithmetic and 
# run the following command:
#
# > ./gunplot4b.sh
#


# Check that the utility is present.

type gnuplot >/dev/null 2>&1 || { echo >&2 "You need gnuplot to run the script. Aborting."; exit 1; }


# Run test programs.

srcdir=../test
$srcdir/etest5b $srcdir/testmat3.mtx rvalues.mtx -e ai -ss 20


# Draw distribution of Ritz values.

filename=rvalues.mtx
while read M N L B X; do
    echo $M | grep -v '^%.*' > /dev/null
    if [ $? -eq 0 ]; then
	size=$M
	break
    fi
done < $filename

gnuplot -e "filename='$filename'; size='$size'" gnuplot4b.plt

rm rvalues.mtx
