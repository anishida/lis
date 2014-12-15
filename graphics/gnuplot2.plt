# Sample script to draw pattern of a square matrix

set terminal postscript eps enhanced color
set output "matrix_pattern.eps"

set title 'Pattern of Matrix'
set pm3d
set pm3d map
unset key
unset xtics
unset ytics
unset colorbox


# Read data in Matrix Market format. 

set datafile commentschars "%"
set size square
set xrange [1:size] 
set yrange [1:size] reverse 


# Draw pattern of nonzero entries.

splot filename every ::1::nnz with points pointtype 7 pointsize 1
