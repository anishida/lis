# Sample script to draw Ritz value distribution on complex plane

set terminal postscript eps enhanced color
set output "Ritz_value_distribution.eps"

set title "Distribution of Ritz Values on Complex Plane"
set xlabel "x"
set ylabel "y"
set style line 1
unset key


# Read data in Matrix Market format. 

set datafile commentschars "%"


# Draw Ritz value distribution.

plot filename using 2:3 with points pointtype 7 pointsize 1
