# Sample script to draw eigenvectors

set terminal postscript eps enhanced color
set output "eigenvectors.eps"

set title "Eigenvectors"
set xlabel "Vector Component"
set ylabel "Component Value"
set style line 1
unset key


# Read data in Matrix Market format. 

set datafile commentschars "%"
set xrange [0:size] 


# Draw eigenvectors.

plot for [i=1:subspace] filename every subspace::i::size*subspace using 1:3 with lines
