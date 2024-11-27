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
# set logscale x
# set logscale y

# plot "obs/absdetPhi1.000000.dat" using 0:2
# replot "obs/absdetPhi2.200000.dat" using 0:2
# replot "obs/absdetPhi3.400000.dat" using 0:2
# replot "obs/absdetPhi4.600000.dat" using 0:2

binwidth=0.1
bin(x,width)=width*floor(x/width)
plot "obs/absdetPhi1.000000.dat" using (bin($2,binwidth)):(1.0) smooth freq with boxes
replot "obs/absdetPhi2.200000.dat" using (bin($2,binwidth)):(1.0) smooth freq with boxes
replot "obs/absdetPhi3.400000.dat" using (bin($2,binwidth)):(1.0) smooth freq with boxes
replot "obs/absdetPhi4.600000.dat" using (bin($2,binwidth)):(1.0) smooth freq with boxes
# replot [0.001:0.1] x**2
# set output