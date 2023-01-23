module Aglomerar

push!(LOAD_PATH, "../Preprocesamiento/")

using PyPlot
using Statistics, StatsBase
using Clustering
using Distances  # Este es necesario para hacer la matriz de distancia, es complemento a Clustering.

using ArraySetTools  # pa aprovechar las cosas que ya tenemos.

export normalizar1, filtraclusterchicos, declustaset,
     dictatabla, plot4Ddiscs, scatterclust



function normalizar1(xx::Array)
    norma=maximum(abs.(xx))
    result=xx./norma
    return result
end



function filtraclusterchicos(puntos::Array, clustree, umbral)
    tabla=hcat(puntos,clustree)
    valores=unique(clustree)
    dd=Dict([(i,count(x->x==i,clustree)) for i in valores])
    predicado(j)=dd[j]>umbral
    result=tabla[map(x->predicado(x), clustree),:]
    (clusterfuck, cual)=findmax(dd)
    println("el cluster mayor tiene ",clusterfuck, " elementos, es el ", cual)
    return result
end

function filtraclusterchicos(puntos::Dict, umbral)
    # no solo aqui ya estan separados por dicctionario, sino tambien
    # en "electrodos", es decir numeros enteros
   result=Dict{Int, Set}()
    for k in keys(puntos)
        if length(puntos[k])>umbral
            result[k]=puntos[k]
        end
    end
    return result
end




function declustaset(puntosyclust::Array; enteros=true)
    #= la funcion convierte el cluster en una tabla
    a un conjunto de enteros, sobre una malla 
    =#
    
    result=Dict{Int,Set}()
    clustnames=unique(puntosyclust[:,end])
    
    for q in clustnames
        subt=puntosyclust[puntosyclust[:,end].==q,:]
        if enteros
        punt=rowstoset(round.(Int64,subt[:,1:end-1]))
        else
            punt=rowstoset(subt[:,1:end-1], enteros=enteros)
        end    
        qindez=round(Int64,q)
        result[qindez]=punt
    end
    return result

end



function dictatabla(dict::Dict; ancho=4, una=true)
    # si una, entonces todo va en una tabla, si no, va en tablitas separadas
    #funcion solo para plotear, no se si luego la vas a usar.
    
    if una 
        aux=zeros(1,ancho+1)
    else
        aux=Dict{Any, Any}()
    end
     println("vamos bien 1 ")
    
    for q in keys(dict)
        
        aux2=zeros(1,ancho)
        
        for renglon in dict[q]
            aux2=vcat(aux2, renglon)
            rintln("vamos bien 2")
        end
        
        println("vamos bien 3")
        
        aux2=aux2[2:end,:]
        
        (l,k)=size(aux2)
        columnaextra=ones(Int64, l).*q
        aux3=hcat(aux2, columnaextra)
        
        println("a ver")
        
        
        if una
            aux=vcat(aux, aux3)
        else
            aux[q]=aux2
        end
    end
    
    if una
        return aux[2:end,:]
    else
        return aux
    end
end


function plot4Ddiscs(datos::Array; escala=0.05)
    figure(figsize=(8,6.5))
    
    xx=datos[:,1]
    yy=datos[:,2]
    gordis=abs.(datos[:,3])
    nn=datos[:,4]
    
    xlim(0.0,65.0)
    ylim(0.0,65.0)
    
    
    scatter(xx,yy, s=escala*gordis, edgecolors="darkmagenta", c=nn, cmap="plasma", lw=0.5)
    colorbar()
    show()
end


function scatterclust(puntos, clustn)

fafa=figure(figsize=(4,4))
#axis("equal")

xlim(0.0,65.0)
ylim(0.0,65.0)
scatter(puntos[:,1],puntos[:,2], s=0.5, c=clustn, cmap="inferno")
colorbar(fraction=0.045)
end


end #module

