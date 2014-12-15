#!/bin/sh

# Sample shell script to draw matrix pattern using gnuplot
# Run the following command:
#
# > ./gunplot2.sh matrix_filename
#


# Check that the utility is present.

type gnuplot >/dev/null 2>&1 || { echo >&2 "You need gnuplot to run the script. Aborting."; exit 1; }


# Draw matrix pattern.

filename=$1
while read M N L B X; do
    echo $M | grep -v '^%.*' > /dev/null
    if [ $? -eq 0 ]; then
	size=$M
	nnz=$L
	break
    fi
done < $filename

gnuplot -e "filename='$filename'; size='$size'; nnz='$nnz'" gnuplot2.plt



