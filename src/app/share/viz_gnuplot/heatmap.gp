set terminal svg size 600,600 font "Arial,10"
set output "@OUTPUT@"
set view map
unset key
set xrange [0:*]
set yrange [0:*]
set size ratio -1
plot "@DATA@" using 1:2:3 with image notitle
