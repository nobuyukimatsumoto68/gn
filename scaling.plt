# set term postscript eps color blacktext enhanced
# set output 'scale.eps'
set xlabel font "Liberation Sans,16"
set ylabel font "Liberation Sans,16"
set tics font "Liberation Sans,20"
set title font "Liberation Sans,22"
# set title sprintf('Wilson-Majorana, %s=%d', '{/Symbol n}', nu) # "t2 t2 rectangular lattice"
# set lmargin 0.
# set rmargin 0.
# set tmargin 0.
# set bmargin 0.
# set yrange [0.0:1.0]
# set zrange [-6.0:6.0]
# set view 74, 206
# set view 80, 340
# set view 85, 170
set xlabel "x"
set ylabel "y"
# set zlabel ""
set logscale x
set logscale y
plot "scaling.dat" using 1:(abs($2))
replot [0.001:0.1] x**2
# set output