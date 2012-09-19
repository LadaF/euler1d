set terminal jpeg
set output "calc1.jpg"
set multiplot
set size 1,0.5
set origin 0, 0

set nokey

set title "Velocity"
set xrange[-1:1]
set yrange[0:1.5]

plot "calc1.dat" using 1:3 smooth csplines

set title "Density"
set yrange[0:1.1]
set size 1,0.5
set origin 0, 0.5
plot "calc1.dat" smooth csplines

