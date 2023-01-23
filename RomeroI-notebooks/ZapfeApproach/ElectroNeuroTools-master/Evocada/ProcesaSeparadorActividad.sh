#!/bin/bash


if [ $# -eq 0 ]
  then
    echo "voy  usar el directorio /media/discogrande/estimulacion_control2/"
    directoriolugar="/media/discogrande/estimulacion_control2/"
elif
    directoriolugar=$1
fi


for i in $(ls $directoriolugar/*.brw);
do
    j=$(basename $i .brw)
    julia SeparaActividadyPromedia.jl $j
    echo "Uno mas listo"
    echo "\n"
 
done
	 
