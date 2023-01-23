#!/bin/bash

##ploteando un chingo con varios archivos

datosfondobasename="LaplacianFilteredSpaceTime-"
numero=$1
dataname=$datosfondobasename$numero".dat"
#usandonegativo='(-$2/4.):(-$1/4):($3/500.)'
#usandopositivo='($2/4.):($1/4):($3/500.)'
#usandomatrix='($1/4.):($2/4):($3)'

usandonegativo='(-$2):(-$1):($3/1500.)'
usandopositivo='($2):($1):($3/1500.)'
usandomatrix='($1):($2):($3)'


outname=`basename $dataname .dat`.png
index="::::$numero"

#segundos=`echo "scale=10; $numero)/7022.0" | bc`
#titulo="Repolarizacion, Centro de Masa, $segundos"
#croac='$titulo'

echo $outname

gnuplot <<EOF

set xr[-0.50:63.5]
set yr[-0.5:63.5]
set cbr[-10:10]
set size ratio -1
set pointsize 0.35
unset key
set term pngcairo size 400,400 fontscale 0.75
segundos=$numero
titulo=sprintf("Centro de Masa %d cuadros", segundos)
set title titulo
set palette defined (0 "#000077", 1 "#0000FF", 2 "white", 3 "#FF0000", 4 "#990000")
#set palette defined (0 "white", 1 "#0000FF", 2 "black")
set out "$outname"
plot "$dataname" using $usandomatrix matrix  w image, "DatosCMNegativo01.dat" usi $usandonegativo every $index w lp ls 7 palette ,"DatosCMPositivo01.dat" usi $usandopositivo every $index w lp ls 7 palette 
 
set out

EOF
