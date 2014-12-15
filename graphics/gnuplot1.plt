# Sample script to draw residual histories

set terminal postscript eps enhanced color
set output "residual_histories.eps"

set title "Residual Histories"
set xlabel "Number of Iterations"
set ylabel "Relative Residual 2-norm"
set style line 1
set logscale y
set format y "10^{%L}"


# Draw residual histories of test programs.

plot \
"ILU0-BiCG" with linespoints, \
"ILU0-CGS" with linespoints, \
"ILU0-BiCGSTAB" with linespoints
