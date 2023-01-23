#=
Esta Rutina de Julia incorpora las definiciones y cositas hechas en los
notebooks de este directorio para rapidamente estimar cuantos canales
están inusables y como se ven
=#

using JLD
namefile=ARGS[1]
println(" Usando $namefile")
namecola=namefile[end-2:end]
println(namecola)
if namecola=="brw"
    using HDF5
elseif namecola=="jld"
    println("chido")
else
    error("No se que tipo de archivo es la terminación $namecola")
end
    
stringgeneral=namefile[1:end-4]

include("SeparaActividadySaturados01.jl")
using SeparaActividadySaturados01
#ya en el modulo estan las funciones necesarias.

if namecola=="brw"
    Datos=AbreyCheca("$namefile")
    factor=Datos["factor"]
    freq=Datos["frecuencia"]/1000 #cuadros por segundo
#la mayoría de los datos están en una lista monstruosa
    DatosCrudosArreglados=reshape(Datos["DatosCrudos"], (4096, Datos["numcuadros"]));
    LFP=FormaMatrizDatosCentrados(DatosCrudosArreglados, factor)
elseif namecola=="jld"
    Datos=load(namefile)
    freq=Datos["freq"]
    if haskey(Datos, "LFPSaturados")
        LFP=Datos["LFPSaturados"]
    else
        LFP=Datos["LFPTotal"]
    end
end
        

println("Esta es la frecuencia de muestreo: ", freq)
tmax=size(LFP,3)
ini=150
fini=tmax
println("Este es el número de cuadros: ", tmax)


PruebaRespuesta=BuscaCanalRespActPot(LFP,freq, 120,-120,-800)
#primero por promedio
Saturados1=BuscaSaturados(LFP,1200,ini,fini)
# luego por ventanas de desviacion
Saturados2=BuscaSaturadosStd(LFP,50,2)

Saturados=union(Saturados1,Saturados2)
setdiff!(PruebaRespuesta,Saturados)

tantossaturados=length(Saturados)
tantosrespuesta=length(PruebaRespuesta)

println("Me encontre ", tantossaturados, " saturados y ", tantosrespuesta, " posibles respuestas")


outname=string(stringgeneral,".jld")
save(outname,
     "LFPTotal", LFP,
     "freq",freq,
     "Canalesrespuesta", PruebaRespuesta,
     "CanalesSaturados", Saturados)
