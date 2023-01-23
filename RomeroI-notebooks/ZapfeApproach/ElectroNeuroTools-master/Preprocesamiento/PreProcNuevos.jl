push!(LOAD_PATH, ".")
push!(LOAD_PATH, "../Isabel/")
using Z_auxiliaresBRW
using PreprocUInt8
using HDF5, PyPlot, Statistics
using ArraySetTools

arxname="/home/karel/BRWFiles/Isabel2019/Cacho_06_control_02.brw"

if occursin("Cacho", arxname)
    arx=h5open(arxname)
    nota=read(arx["notacacho"])
    print(nota)
else
    arx=brw_things(arxname)
    print("Estas partiendo de un brw original")
end

freq=read(arx["SamplingRate"])/1000 # usamos kHz y milisegundos aqui porfa
dataraw=read(arx["dset"])
factor=read(arx["Factor"])
listchans=read(arx["Chs"])

evocado=false


tamanos=size(dataraw)
if length(tamanos)==2
    ncuadros=tamanos[2]
else
    tot=tamanos[1]
    ncuadros=Int(tamanos/4096)
end

#Estas medidas estan en cuadros, no en ms.
if evocado
    retrazo=round(Int, ceil(5*freq))
    final=round(Int, ceil(60*freq))
    latencia=round(Int, ceil(1.5*freq))
else
    retrazo, final, latencia = 0, 0, 0
end

tiempototalms=round(ncuadros/freq; digits=1) 
println("Tienes ", ncuadros, " cuadros de muestreo a ", round(freq; digits=4), " cuadros por milisegundo")
println( "Esto corresponde a  ", tiempototalms, "ms." )

#los parametros son datos promediados, frecuencia, tiempo post estimulo en ms, umbral en microvolts, umbral de saturación en microvolts

iniciobusqueda=1
finbusqueda=282

PruebaRespuesta=buscaCanalPicos(dataraw, factor=factor,
    tini=iniciobusqueda, tfin=finbusqueda, maxvolt=-100, freq=freq)
# y luego los saturados


satu=buscasatura(dataraw, factor=factor,
                 tini=iniciobusqueda, tfin=finbusqueda,
                 umbral=2000)


numerosaturados=length(satu)
println("Encontramos ", numerosaturados, " canales probablemente saturados.")
println("Encontramos ", length(PruebaRespuesta), " canales probablemente con actividad.")




quietos=buscastdraras(dataraw,
                      tini=iniciobusqueda,tfin=finbusqueda,
                      factor=factor, freq=freq,
                      bajo=16,alto=10000)
intensos=buscastdraras(dataraw, tini=iniciobusqueda, tfin=2*finbusqueda,
                       factor=factor,
                       freq=freq, bajo=0, alto=30)

numerosaturados=length(quietos)
println("Encontramos ", numerosaturados, " canales muy quietos")
numerosaturados=length(intensos)
println("Encontramos ", numerosaturados, " canales muy intensos")


# Inspeccion visual

DesviacionPorCanal=zeros(64,64)
for j=1:64
    for k=1:64
        ChorizoExemplo=vec(dataraw[k+(j-1)*64,:])
        DesviacionPorCanal[j,k]=std(ChorizoExemplo)
    end
end


ion()
fff=figure(figsize=(8,8))
xlim(0,65)
ylim(0,65)
title("σ")
limites=30
imagen=imshow(DesviacionPorCanal, origin="lower", interpolation="nearest",cmap="gnuplot2", 
                      vmin=0,vmax=limites, extent=[0.5,64.5,0.5,64.5])
cb=colorbar()

x=[]
y=[]
for j in PruebaRespuesta
  append!(x,j[3])
   append!(y,j[2])
end


x2=[]
y2=[]
for j in satu
  append!(x2,j[3])
   append!(y2,j[2])
end


x3=[]
y3=[]
for j in quietos
  append!(x3,j[3])
   append!(y3,j[2])
end


x4=[]
y4=[]
for j in intensos
  append!(x4,j[3])
   append!(y4,j[2])
end

#=
x5=[]
y5=[]
for j in ruidosos
  append!(x5,j[2])
   append!(y5,j[1])
end
=#
scatter(x,y, marker="o",c="lightblue", s=15)
scatter(x2,y2, marker="o",c="green", s=25)
scatter(x3,y3, marker="x",c="green", s=25)
scatter(x4,y4, marker="+",c="red", s=25)
#scatter(x5,y5, marker="v",c="green", s=25)

println( "el tipo de intensos ", typeof(intensos))
println(" el tipo de PruebaRespuesta ", typeof(PruebaRespuesta))
mal=union(quietos, satu)
buenos=union(intensos, PruebaRespuesta)
setdiff!(buenos, mal)

println(" el tipo de buenos ", typeof(buenos))

malforsave=elemtorow(mal, n=3)
buenforsave=elemtorow(buenos, n=3);

println(" el tipo de buenforsave ", typeof(buenforsave))

dirgen=dirname(arxname)*"/"
basegen=basename(arxname)[1:end-4]
#una palabra para indicar el intervalo 
outname=string(dirgen, basegen, "_preproc.h5")
#aqui hay que arreglar el rollo del directorio

listaextra=Dict(
     "CanalesBuenos" => buenforsave,
     "CanalesMalos" => malforsave)

println(outname )

   h5open(outname, "w")  do file
        for nam in names(arx)
            datos=read(arx[nam])
            println("abriendo los datos ", nam)
    
            write(file, "$nam", datos)
        end
        
        for k in keys(listaextra)
            datos=listaextra[k]
            println("abriendo los datos ", k)
            println
            write(file, "$k", datos)
        end
    end

println("hemos guardado este cacho en el archivo ", outname)


readline()
