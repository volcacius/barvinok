set output 'sum_p.eps'
set terminal postscript eps 18
set key bottom
set logscale y
set ylabel "time (s)"
plot "sum_p_times_0" title "reduction to unweighted Barvinok", "sum_p_times_2" title "Euler-Maclaurin", "sum_p_times_3" title "Laurent expansion (old)", "sum_p_times_4" title "Laurent expansion (new)"
