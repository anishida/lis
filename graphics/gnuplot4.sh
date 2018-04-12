#!/bin/sh

# Sample shell script to draw eigenvalue distribution on complex plane
# using gnuplot
# Make executable files enabling complex arithmetic and 
# run the following command:
#
# > ./gunplot4.sh
#


# Check that the utility is present.

type gnuplot >/dev/null 2>&1 || { echo >&2 "You need gnuplot to run the script. Aborting."; exit 1; }


# Run test programs.

srcdir=../test
$srcdir/etest5 $srcdir/testmat3.mtx evalues.mtx /dev/null /dev/null /dev/null -e ai -ss 20


# Draw eigenvalue distribution.

filename=evalues.mtx
while read M N L B X; do
    echo $M | grep -v '^%.*' > /dev/null
    if [ $? -eq 0 ]; then
	size=$M
	break
    fi
done < $filename

gnuplot -e "filename='$filename'; size='$size'" gnuplot4.plt

rm evalues.mtx
