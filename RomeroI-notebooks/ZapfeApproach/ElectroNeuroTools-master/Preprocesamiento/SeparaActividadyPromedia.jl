#=
Programa-script que separa la actividad alrededor de 3 trancazos, 
la centra y la promedia y la guarda en un archivo jld. Ademas,
selecciona los probables canales saturados y los probables canales
con respuesta fisiológica, y los guarda en el mismo archivo como
conjunos de datos
=#

include("SeparaActividadySaturados01.jl")

using SeparaActividadySaturados01
using JLD


if length(ARGS)==0
    println("Dar el nombre de un archivo para trabajar")
else
    abrestring=ARGS[1]
    stringgeneral=replace(abrestring, ".brw", "")
    println("Voy a trabajar con el archivo ", abrestring)
end
    

#stringgeneral="est2_con3"

#Todos los datos de la matriz del brw file
Datos=AbreyCheca(abrestring)

#Los ponemos bien, array espacio-temporal
if size(Datos["DatosCrudos"])[1] != 4096
    DatosCrudosArreglados=reshape(Datos["DatosCrudos"], (4096, Datos["numcuadros"]))
else
DatosCrudosArreglados=Datos["DatosCrudos"]
end

#parametros globales.
freq=Datos["frecuencia"]/1000 #cuadros por milisegundo.
factor=Datos["factor"] #Factor de conversion de numeros enteros a microVolts

retrazo=round(Int, ceil(5*freq))
final=round(Int, ceil(60*freq))
latencia=round(Int, ceil(1.5*freq))
#Gis recomienda esperar hasta 6.5 ms por cada estimulo.
tiempopostgolpe=round(Int,ceil(6.5*freq))

exemplo=vec(DatosCrudosArreglados[1500,:]);
listongas=EncuentraTrancazosRaw(exemplo,500)
println("Sitios de trancazo electrico=", listongas)

ActividadRaw=ActivAlrededorTrancazo(listongas, DatosCrudosArreglados, freq)
(bla,tmax)=size(ActividadRaw["Trancazo_1"])
println("Tiempo, en cuadros, de los intervalos con trancazo=",tmax)


DatosCentrados=Dict{AbstractString, Array}()

for k in keys(ActividadRaw)
    DatosCentrados[k]=FormaMatrizDatosCentrados(ActividadRaw[k], factor)
end

LFPPromedio=(DatosCentrados["Trancazo_1"]+DatosCentrados["Trancazo_2"]+DatosCentrados["Trancazo_3"])/3


PruebaRespuesta=BuscaCanalRespActPot(LFPPromedio,freq, 120,-120,-800)

iniciobusqueda=retrazo
finbusqueda=retrazo+4*latencia
Saturados=BuscaSaturados(LFPPromedio,1200,iniciobusqueda,finbusqueda)
numerosaturados=length(Saturados)

println("Encontramos ", numerosaturados, " canales probablemente saturados.")


outname=string("LFP_Promedio_",stringgeneral,".jld")
parametros=(freq, retrazo, latencia,final)
save(outname,
     "LFPSaturados", LFPPromedio,
     "freq",freq,
     "retrazo", retrazo,
     "latencia", latencia,
     "final", final,
     "Canalesrespuesta", PruebaRespuesta,
     "CanalesSaturados", Saturados)


