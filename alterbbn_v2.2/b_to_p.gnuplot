set datafile separator ','

set key autotitle columnhead

set ylabel "Baryon-to-photon ratio"
set xlabel "Plasma Temperature (MeV)"
set logscale x 10
set xrange[15:0.0007]
set logscale y 10
#set yrange[1e-15:1e-7]

set xtics nomirror
set ytics nomirror
set key right top

set style line 1 lc rgb "#e6d32e" lw 2 #gold for std cosmo model
set style line 2 lc rgb "#782ee6" lw 2 #purple for our model
set style line 3 lc rgb "#e6872e" lw 2 #orange for phi model

plot 'evolution.csv' using 3:21 with lines ls 1, 'evolution_vs.csv' using 3:21 with lines ls 2, 'evolution_phi.csv' using 3:21 with lines ls 3
