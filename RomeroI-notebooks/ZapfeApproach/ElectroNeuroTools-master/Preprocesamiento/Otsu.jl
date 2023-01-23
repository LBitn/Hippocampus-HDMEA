module Otsu

using Statistics
using StatsBase

export stdventana, otsumethod, otsuumbralizar




function stdventana(datos::Array, ancho::Int, paso::Int)
    # de un array lineal, sacar la desviacion
    # estandar por ventanas de ancho dado
    anchomedio=floor(ancho/2)
    tantos=floor(Int, (length(datos)-ancho)/paso)
    result=zeros(tantos)
    for t=1:tantos
        result[t]=std(datos[(t-1)*paso+1:(t-1)*paso+ancho])
    end
    return result
end

function otsumethod(datos::Array)
    #=
    Metodo de Otsu para separar una imagen en dos subconjuntos.
    =#
    minimo=minimum(datos)
    maximo=maximum(datos)
    binsdefault=2*ceil(Int,sqrt(length(datos)))
    h=fit(Histogram, vec(datos),closed=:right,nbins=binsdefault)
    
    (rango, cuentas)=(h.edges[1], h.weights) #stupid Statistics.
    tantos=length(rango)
    #valores
    omega1=0
    omega2=0
    mu1=0
    mu2=0
    sigmab=0
    sigmabtemp=0
    tbest=0
    varlim=0
    for t=1:tantos-1
        omega1=sum(cuentas[1:t])
        omega2=sum(cuentas[t+1:tantos-1])
        mu1=sum(cuentas[1:t].*rango[1:t])/omega1
        mu2=sum(cuentas[t+1:tantos-1].*rango[t+1:tantos-1])/omega2        
        sigmabtemp=omega1*omega2*((mu1-mu2)^2)
        if sigmabtemp>sigmab
            sigmab=sigmabtemp
            tbest=t
            varlim=rango[t]
        end
    end
    return (sigmab,tbest,varlim)
end


function otsuumbralizar(datos::Array)
    # aplanar datos
    result=zeros(size(datos))
    umbral=otsumethod(datos)[3]
    result=map(x->(x>umbral) ? 1 : 0, datos)
    return result
end


end

