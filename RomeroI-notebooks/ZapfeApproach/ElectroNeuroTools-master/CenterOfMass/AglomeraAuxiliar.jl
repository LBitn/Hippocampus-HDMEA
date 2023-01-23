""" Modulo que contiene la talacha rutinaria del aglomeramiento """

module AglomeraAuxiliar

push!(LOAD_PATH, ".")

using Aglomerar
using StatsBase # para nquantil
using Distances # para pairwise euclidean matrix


export dictatablaord, filtraquantil, hazmatrizdist

""" funcion de dict a tabla ordenada """
function dictatablaord(dict::Dict)
    
    fruncio=Array{Float64}(undef, 0,4)

    for k in keys(dict)
        dat=dict[k]
        (alto, ancho)=size(dat)
        auxt=ones(alto).*k
        afafa=hcat(dat, auxt)
        fruncio=vcat(fruncio,afafa)
    end

    fruncio=fruncio[sortperm(fruncio[:, 4]), :];

    return fruncio

end

function filtraquantil(datos::Array; cual=3)
    # gmin:= gordo minimo
    gmin=nquantile(datos[:,cual], 10)[2]
    result=datos[datos[:,cual] .> gmin, :]
    return result
end


function hazmatrizdist(datos::Array)
    datosaux=permutedims(datos)
    datosaux[3,:]=normalizar1(datosaux[3,:]).*64
    datosaux[4,:]=normalizar1(datosaux[4,:]).*64
    result=pairwise(Euclidean(), datosaux, dims=2)
    return result

end
    
end #module
