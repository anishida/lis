# Sample script to draw distribution of Ritz values on complex plane

set terminal postscript eps enhanced color
set output "distribution_of_Ritz_values.eps"

set title "Distribution of Ritz Values"
set xlabel "Real part"
set ylabel "Imaginary part"
set style line 1
unset key


# Read data in Matrix Market format. 

set datafile commentschars "%"


# Draw distribution of Ritz values.

plot filename using 2:3 with points pointtype 7 pointsize 1
