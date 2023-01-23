#= funciones para ayudar a trabajar con las trayectorias
Contiene funciones que reducen el catenario (el diccionario
de trayectorias) y que quitan puntos "antes que" o "despues que"
cierto tiempo
=#

module TrayectoriasAux

using Statistics

export reducecatenario, AntesQue, DespuesQue, distprom,
    dist, set2dict, cortacatenarios

dist(x1,y1,x2,y2)=sqrt((x1-x2)^2+(y1-y2)^2)

function set2dict(xx::Set)
    j=1
    result=Dict{Integer, Any}()
    for x in xx
        result[j]=x
        j+=1
    end
    return result
end

function reducecatenario(Catenario::Dict;
                         gordmin=25,minlong=7,maxlong=100, distmin=0)
    
    result=Dict{Integer, Array{Any}}()
    
    for (k,p) in Catenario
    
        distrec=dist(p[1,1],p[1,2],p[end,1],p[end,2])
        
        gordura=abs.(p[:,3])
        longus,gordus=size(p)
        if (mean(gordura)>gordmin) &&
            (longus>minlong) &&
            longus < maxlong && distrec>distmin
            result[k]=p
        end
    end

    return result
end


function cortacatenarios(Catenario::Dict; nini,nfin)
    # reduce las trayectorias por cuadro de inicio y final

    result=Dict{Integer, Array{Any}}()
    
    for (k,p) in Catenario

        paux=AntesQue(p,nfin)
        paux=DespuesQue(paux,nini)
        if length(paux) != 0
            result[k]=paux
        end
         
   end

    return result
end


function AntesQue(Datos::Array, tiempo)
    Cadena=[0 0 0 0]
 
    for a in eachindex(Datos[:,4])
        if Datos[a,4]<tiempo
  
            Cadena=vcat(Cadena, transpose(Datos[a,:]))
          
    end
end
    Cadena=Cadena[2:end,:]
    return Cadena
end


function DespuesQue(Datos::Array, tiempo)
    Cadena=[0 0 0 0]
    for a in eachindex(Datos[:,4])
        if Datos[a,4]>tiempo
            Cadena=vcat(Cadena, transpose(Datos[a,:]))
    end
end
    Cadena=Cadena[2:end,:]
    return Cadena
end

function distprom(Catenario::Dict)
    # Saca la velocidad promedio por trayectoria
    result=Dict{Integer, Float64}()
    for (k,p) in Catenario
        longus,anchos=size(p)
        if longus>1
            d=0
            for j=1:longus-1
                d+=dist(p[j,1],p[j,2],p[j+1,1],p[j+1,2])
            end
            d/=(longus-1)
        end
        result[k]=d
    end
    return result
end

function distprom(Catenario::Set)
    aux=set2dict(Catenario)
    distprom(aux)
end

end # module
