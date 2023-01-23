#=
Un standalone solo para probar parametros de busqueda de saturados.
=#

using JLD
include("SeparaActividadySaturados01.jl")
using SeparaActividadySaturados01

if length(ARGS)==0
    error("Dar el nombre de un archivo LFP_*.jld para graficar")
else
    archivo=ARGS[1]
    println("Voy a trabajar con el archivo ", archivo)
end
    

#Los datos para fraficar
LFP=load(archivo)["LFPSaturados"]
saturados=load(archivo)["CanalesSaturados"]
respuestas=load(archivo)["Canalesrespuesta"];
frecuencia=load(archivo)["freq"]
retrazo=load(archivo)["retrazo"]


palabra=replace(archivo,".jld", "")
println(palabra)
palabrita=replace(palabra,"otroIntento", "")

if !(isdir(palabrita))
     mkdir(palabrita)
end

latencia=round(Int, ceil(1.5*frecuencia))

iniciobusqueda=retrazo
finbusqueda=retrazo+10*latencia

println("Estamos buscando desde el cuadro ",iniciobusqueda, " hasta el ",
        finbusqueda)

saturados=BuscaSaturados(LFP,1100,iniciobusqueda,finbusqueda)

numerosaturados=length(saturados)

println("Encontramos ", numerosaturados, " canales probablemente saturados.")

tutu=load(archivo)
tutu["CanalesSaturados"]=saturados
save(archivo, tutu)




