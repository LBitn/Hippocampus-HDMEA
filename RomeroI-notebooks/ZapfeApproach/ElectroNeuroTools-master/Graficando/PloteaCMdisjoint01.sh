#!/bin/bash

##ploteando un chingo con varios archivos

datosfondobasename="CSD-"
datosCMNegname="CMNegLapla-"
datosCMPosname="CMPosLapla-"
numero=$1
dataname=$datosfondobasename$numero".dat"
datosCMNeg=$datosCMNegname$numero".dat"
datosCMPos=$datosCMPosname$numero".dat"
#usandonegativo='(-$2/4.):(-$1/4):($3/1000.)'
#usandopositivo='($2/4.):($1/4):($3/1000.)'
#usandomatrix='($1/4.):($2/4):($3)'

usandonegativo='($1):($2):(abs($3*0.001))'
usandopositivo='($1):($2):($3*0.001)'
usandomatrix='($1):($2):($3)'

outname=`basename $dataname .dat`ConFondo.png


#segundos=`echo "scale=10; $numero)/7022.0" | bc`
#titulo="Repolarizacion, Centro de Masa, $segundos"
#croac='$titulo'

echo $outname

gnuplot <<EOF

set xr[-0.50:63.5]
set yr[63.5:-0.5]
set cbr[-150:150]
set size ratio -1
set pointsize 0.35
unset key
set term pngcairo size 600,600 fontscale 1.25 font "Inconsolata"
segundos=$numero/7022.0
titulo=sprintf("Centre of Mass,: %3.5f s ", segundos)
set title titulo
set palette defined (0 "#000077", 1 "#0000FF", 2 "white", 3 "#FF0000", 4 "#990000")
unset cbtics
set label "Source" at graph 1.1, 0.98
set label "Sink" at graph 1.1, 0.02
set style fill transparent solid 0.75 noborder
#set palette defined (0 "white", 1 "#0000FF", 2 "black")
set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')
set out "$outname"
#set pm3d interpolate 4,4 map explicit
set dgrid 256,256 splines
set table "Tabla.dat"
splot "$dataname" using $usandomatrix matrix 
unset table
set pointsize 1
#unset pm3d

plot "Tabla.dat" using $usandomatrix   w image, \
"$datosCMNeg" usi $usandonegativo  w circles lw 3 lc rgb "#4a90ff",\
"$datosCMPos" usi $usandopositivo  w circles lw 3 lc rgb "#ff1414",\
"$datosCMNeg" usi $usandonegativo  w points ls 7 lc rgb "#0000ff",\
"$datosCMPos" usi $usandopositivo  w points ls 7 lc rgb "#ff0000"


#plot  "$datosCMNeg" usi $usandonegativo  w circles lw 3 lc rgb "#000075", "$datosCMPos" usi $usandopositivo  w circles lw 3 lc rgb "#750000"
set out

EOF
