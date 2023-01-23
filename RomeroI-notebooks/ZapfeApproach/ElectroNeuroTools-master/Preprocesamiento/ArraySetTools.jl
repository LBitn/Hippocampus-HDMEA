module ArraySetTools

export rowstoset, elemtorow, orderelectrodes


## ******************************************
## funciones para manipular arrays y sets



function rowstoset(xxs::Array; enteros=true)
    # pasa los renglones de un arreglo a un conjunto como si jueran elementos
    if enteros
        result=Set{Array{Int64,2}}()
    else
        result=Set{Array{Float64,2}}()
    end
    (nrow,ncol)=size(xxs)
    for j in 1:nrow
        push!(result, reshape(xxs[j,:], 1, ncol))
    end
    return result
end

function elemtorow(xxs::Set; n=2, renglones=false)
    #pasa los elementos de un conjunto a rengoles de un array
    result=zeros(UInt16, 1,n)
    for k in xxs
        if renglones
            result=vcat(result, k)
        else
            result=vcat(result, transpose(k))
         end
    end
    return result[2:end,:]
end


#=
function puntoZ2enlinea(puntoini::Array, direccion)
    # los puntos estan en "renglon columna"
    xo=puntoini[2]
    yo=puntoini[1]
    
    
end

function orderelectrodes(xxs::Set, puntoinicial::Array;
                         dextro=true, inidir=3.1415)
    

end

=#


end #module
