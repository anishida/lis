#!/bin/sh

# Sample shell script to draw eigenvalue distribution using gnuplot
# Make executable files and run the following command:
#
# > ./gunplot3.sh
#


# Check that the utility is present.

type gnuplot >/dev/null 2>&1 || { echo >&2 "You need gnuplot to run the script. Aborting."; exit 1; }


# Run test programs.

srcdir=../test
$srcdir/etest5 $srcdir/testmat.mtx evalues.mtx /dev/null /dev/null /dev/null -e si -ie ii -ss 100


# Draw eigenvalue distribution.

filename=evalues.mtx
while read M N L B X; do
    echo $M | grep -v '^%.*' > /dev/null
    if [ $? -eq 0 ]; then
	size=$M
	break
    fi
done < $filename

gnuplot -e "filename='$filename'; size='$size'" gnuplot3.plt

rm evalues.mtx
