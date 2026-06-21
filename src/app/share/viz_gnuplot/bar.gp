set terminal svg size 760,420 font "Arial,11"
set output "@OUTPUT@"
set title "@CAPTION@"
set style data histograms
set style fill solid 0.82 border lc rgb "#2b2b2b"
set boxwidth 0.62
set grid ytics lc rgb "#d9e2e8" lw 1
set key off
set xlabel "@XLABEL@"
set ylabel "@YLABEL@"
set border lc rgb "#4b5563"
set tics textcolor rgb "#374151"
set title font "Arial,14"
set xlabel font "Arial,11"
set ylabel font "Arial,11"
set xzeroaxis lw 1.5 lc rgb "#111827"
set xtics rotate by -25
plot "@DATA@" using 2:xtic(1) lc rgb "#2f6f9f" notitle
