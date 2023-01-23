set term cairolatex size 15cm,12cm standalone lw 3
set out "IntentoInset01.tex"
set multiplot
set border 3 lw 3
set xl 'time [ms]'
set yl 'recorded potential [$\mu$V]'
set xr[0:0.145]
set sty dat lines
set tics nomirror
plot "EjemploPlaticaCrudo01.dat" usi ($0/7022.):1 t "raw data", "EjemploPlaticaSuave01.dat"  usi ($0/7022.):1 t "time smoothed signal" lw 2 lc rgb "blue"
set origin 0.16, 0.12
set size 0.3,0.3
set border 15
set tics scale 1
unset key
unset xl
unset yl
set xr[0.04:0.045]
set xtics 0.04,0.0025
set ytics -500,50,300
plot "EjemploPlaticaCrudo01.dat" usi ($0/7022.):1 t "raw data", "EjemploPlaticaSuave01.dat"  usi ($0/7022.):1 t "time smoothed signal" lw 2 lc rgb "blue"
unset multiplot
set out

