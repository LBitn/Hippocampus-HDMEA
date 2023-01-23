#!/bin/bash

##ploteando un chingo

dataname=$1
outname=`basename $1 .dat`.png
unmbral=30
usando='using 1:2:(criterio($3))'

echo $outname

gnuplot <<EOF

set xr[0:65]
set yr[0:65]
set cbr[-150:150]
set size ratio -1

set term pngcairo size 300,300 fontscale 0.5
unset key
set palette defined (0 "#000077", 1 "#0000FF", 2 "white", 3 "#FF0000", 4 "#990000")
set out "$outname"
criterio(x)=(abs(x)<60)?0:x
plot "$dataname" matrix $usando w image
set out

EOF
