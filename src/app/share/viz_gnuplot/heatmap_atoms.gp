set terminal @FORMAT@ size 700,600
set output "@OUTPUT@"
set view map
unset key
set size ratio -1

set bmargin at screen 0.1
set lmargin at screen 0.1
set rmargin at screen 0.95
set tmargin at screen 0.95

@GRID@
@LABELS@

# Границы сетки по данным (колонки 1 и 2 — это x и y)
stats "@DATA@" using 1:2 nooutput

xmin = STATS_min_x
xmax = STATS_max_x
ymin = STATS_min_y
ymax = STATS_max_y

# Если x,y — центры ячеек, то правильный охват: ±0.5
# Диапазон выводится по boxxyerrorbars из данных.

# Статистика данных из файла (третья колонка)
stats "@DATA@" using 3 nooutput

DATA_min = STATS_min
DATA_max = STATS_max

set cbrange [DATA_min : DATA_max]

# Автоматическая генерация тиков
intervals = 5
step = (DATA_max - DATA_min) / intervals

unset cbtics
set cbtics ("0" 0)

do for [i=0:intervals] {
    x = DATA_min + i * step
    if (abs(x) > 1e-9) {
        eval sprintf('set cbtics add ("%.3g" %g)', x, x)
    }
}

# Рисуем
plot "@DATA@" using 1:2:(0.5):(0.5):3 with boxxyerrorbars fs solid 0.1 border lc rgb "#444444" lw 0.1 palette notitle
