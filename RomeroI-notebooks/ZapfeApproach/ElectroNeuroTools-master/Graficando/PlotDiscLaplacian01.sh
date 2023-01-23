#!/bin/bash

##ploteando un chingo

dataname=$1
outname=`basename $1 .dat`.png
#umbral=30
#usando='using 1:2:(criterio($3))'

echo $outname

gnuplot <<EOF

set xr[0.5:62.5]
set yr [0.5:62.5]
set cbr[-500:500]
set size ratio -1

set term pngcairo size 350,350 fontscale 0.75 font "Comic Sans MS"
unset key
set palette defined (0 "#000077", 1 "#0000FF", 2 "white", 3 "#FF0000", 4 "#770000")
#set palette defined (0 "white", 1 "#0000FF", 2 "black")
set out "$outname"
plot "$dataname" matrix  w image
set out

EOF
