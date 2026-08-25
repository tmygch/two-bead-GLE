set terminal pngcairo size 900,700 enhanced
set output "toy6_mu_t.png"

set size 0.62,0.62
unset grid
set logscale

set label 1 at graph -0.240, 0.95 "(a)"
set tics scale 0.7,0.3

set xlabel "time"
set ylabel "memory kernel"

set key at graph 0.35,0.60
set key spacing 1.3 width 0.0
set key samplen 2

set format x "10^{%L}"
set format y "10^{%L}"

plot [0.009:10][0.001:2] \
"N6_TAG0_5/et_vcor.dat" u 1:4 pt 4 ps 0.5 lc rgb "black" ti "sim (1,6)", \
"N6_TAG1_5/et_vcor.dat" u 1:4 pt 6 ps 0.6 lc rgb "red"   ti "sim (2,6)", \
"N6_TAG4_5/et_vcor.dat" u 1:4 pt 8 ps 0.6 lc rgb "blue"  ti "sim (5,6)", \
"Toy6_N6_mu(t)_i001_j006_genth.dat" w l lw 5 lc rgb "black" ti "theory (1,6)", \
"Toy6_N6_mu(t)_i002_j006_genth.dat" w l lw 5 lc rgb "red"   ti "theory (2,6)", \
"Toy6_N6_mu(t)_i005_j006_genth.dat" w l lw 5 lc rgb "blue"  ti "theory (5,6)"