# Sample script to draw eigenvalue distribution on complex plane

set terminal postscript eps enhanced color
set output "distribution_of_complex_eigenvalues.eps"

set title "Distribution of Complex Eigenvalues"
set xlabel "Real part"
set ylabel "Imaginary part"
set style line 1
unset key


# Read data in Matrix Market format. 

set datafile commentschars "%"


# Draw eigenvalue distribution.

plot filename using 2:3 with points pointtype 7 pointsize 1
