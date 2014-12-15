# Sample script to draw eigenvalue distribution

set terminal postscript eps enhanced color
set output "eigenvalue_distribution.eps"

set title "Distribution of Eigenvalues"
set xlabel "Number of Eigenvalues"
set ylabel "Eigenvalue"
set style line 1
unset key


# Read data in Matrix Market format. 

set datafile commentschars "%"
set xrange [0:size] 


# Draw eigenvalue distribution.

plot filename every ::1::size with points pointtype 7 pointsize 1
