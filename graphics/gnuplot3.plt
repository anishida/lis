# Sample script to draw distribution of real eigenvalues

set terminal postscript eps enhanced color
set output "distribution_of_real_eigenvalues.eps"

set title "Distribution of Real Eigenvalues"
set xlabel "Number of Eigenvalues"
set ylabel "Eigenvalue"
set style line 1
unset key


# Read data in Matrix Market format. 

set datafile commentschars "%"
set xrange [0:size] 


# Draw distribution of real eigenvalues.

plot filename every ::1::size with points pointtype 7 pointsize 1
