set output 'sum.eps'
set terminal postscript eps 18
set key bottom
set logscale y
set ylabel "time (s)"
plot "sum_times_0" title "reduction to unweighted Barvinok", "sum_times_1" title "nested sums", "sum_times_2" title "Euler-Maclaurin", "sum_times_3" title "Laurent expansion (old)", "sum_times_4" title "Laurent expansion (new)"
