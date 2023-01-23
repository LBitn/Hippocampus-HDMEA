module ActividadSutil

push!(LOAD_PATH, ".")

using Otsu
using PreprocTools
using Statistics

#=
Modulo para reconocer canales activos o saturados mas sutiles,
e.g. rebanadas de raton
=#


export buscasigmaactiva

deffreq=17.85550205219098 # frecuencia por default en kHz.

function buscasigmaactiva(datos::Array; malos=Set(), venms=50, freq=deffreq)
    # busca canales con una estadistica diferente del resto.
    # usa criterio de Otsu sobre el maximo de las desv. estandar por ventanas.111
    (alto, ancho, largo)=size(datos)
    pur=zeros(alto,ancho)
    for j=1:alto, k=1:ancho
        aux=desviacionventanas(datos[j,k,:], venms, freq)
        pur[j,k]=maximum(aux)
    end
    med=median(pur)
    for q in malos
        reng,col=q[1], q[2]
        pur[reng,col]=med
    end
    return pur
end

end



