#!/bin/bash

##ploteando un chingo

dataname=$1
outname=`basename $1 .dat`.png
#umbral=30
usando='using ($1/4.0)):($2/4):(-$3))'

echo $outname

gnuplot <<EOF

set xr[0.5:63.5]
set yr[0.5:63.5]
set cbr[-2:2]
set size ratio -1

set term pngcairo size 400,400 fontscale 1.0
unset key
unset tics
set label 1 at screen 0.8,0.9 "Fuente"
set label 2 at screen 0.8,0.1 "Pozo"
#set palette rgbform 33,13,10
#set palette defined (0 "#000077", 1 "#0000FF", 2 "white", 3 "#FF0000", 4 "#770000")
#set palette defined (0 "white", 1 "#0000FF", 2 "black")
set out "$outname"
plot "$dataname" matrix  w image
set out

EOF
