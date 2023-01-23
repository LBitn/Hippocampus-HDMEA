#= Programa para convertir  los datos con estrucuras
"cool" a estructuras convencionales en hdf5 y volverlo a guardar.
Intentemos programarlo con Julia 1.0.0
=#

using JLD, HDF5

numarg=length(ARGS)
if numarg<1
    error("Dame un nombre de archivo JLD que tenga LFP guardados")
else
    nombre=ARGS[1] 
end



nombre=ARGS[1]
nom=nombre[1:end-3]
nomout=nom*"h5"


arx=load(nombre)
listanombres=keys(arx)



function deSetaArray!(xname::String, datostotal::Dict)
    println("Tenemos boey un dato latoso de los de " , xname)
    datos=datostotal[xname]
    aux1=Array[[0::Int64,0::Int64]];
    for i in datos
        push!(aux1, i)
    end
    aux2=permutedims(hcat(aux1[2:end,:]...), [2,1])
    aux3=convert(Array{Int64, 2 }, aux2)
    datostotal[xname]=aux3
    return "¡Has modificado el dato catalogado como "*xname*"!"
    
end


function deDictaArray!(xname::String, datostotal::Dict)
    println("Tenemos boey un dato latoso de los de " , xname)
    datos=datostotal[xname]
    aux1=zeros(Float64, (1,4))
    lista=collect(keys(datos))
    sort!(lista)
    
    for k in lista
        aux=datos[k]
        colextra=ones(Float64, size(aux)[1]).*k
        aux2=hcat(colextra,aux)
        aux1=vcat(aux1, aux2)
    end

    aux1=aux1[2:end,:]
    datostotal[xname]=aux1
    return "¡Has modificado el dato catalogado como "*xname*"!"

end
    
    
if in("Canalesrespuesta", listanombres)
    deSetaArray!("Canalesrespuesta", arx)
end

if in("CanalesSaturados", listanombres)
    deSetaArray!("CanalesSaturados", arx)
end

if in("CMPNeg", listanombres)
    deDictaArray!("CMPNeg", arx)
end

if in("CMNeg", listanombres)
    deDictaArray!("CMNeg", arx)
end

if in("CMPos", listanombres)
    deDictaArray!("CMPos", arx)
end

if (haskey(arx, "Canalesrespuesta") ||  haskey(arx, "CanalesSaturados"))
   h5open(nomout, "w")  do file
        for k in keys(arx)
            datos=arx[k]
            println("abriendo los datos ", k)
            println
            write(file, "$k", datos)
        end
    
        #close(fid)   
    end
else
   println("No hay ninguno de los datos en conjunto conflictivos en este archivo.")   
end
