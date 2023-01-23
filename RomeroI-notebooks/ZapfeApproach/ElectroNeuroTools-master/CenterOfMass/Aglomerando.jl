"""
Programa que hace el aglomeramiento jerarquico de Centros de Masa

"""



push!(LOAD_PATH, ".")

using JLD
using AglomeraAuxiliar, Aglomerar
using Clustering

numarg=length(ARGS)

helpdocstring=" Uso: al menos se requiere un argumento: \n
- un nombre de archivo JLD que tenga CM.
Argumentos optativos segundo y tercero:
- h:= el corte del arbol jerarquico en unidades de los electrodos (distancia)\n
- filt:= el numero de electrodos minimo que tiene que
cubrir un aglomerado para ser tomado en cuenta.\n"

if numarg<1
    error(helpdocstring)
else
    nombre=ARGS[1] 
end

nombregeneral=nombre[1:end-4]

nombresalida=string(nombregeneral, "-Aglo.jld")

println("este es el nombre salida: ", nombresalida)


if numarg==1
  evocada=false
    h=1.05
    filt=3
elseif numarg==2
    evocada=ARGS[2]
    h=1.05
    filt=3
elseif numarg==3
    evocada=ARGS[2]
    h=ARGS[3]
    filt=3
elseif numarg==4
    h=ARGS[3]
    filt=ARGS[4]
end

nota= " los parametros son: \n
evocada: $evocada  \n
corte del arbol h:  $h \n
filtro de aglomerado en electrodos filt: $filt"

println(nota)

archivo=load(nombre)
DatosCMP=archivo["CMP"]
DatosCMN=archivo["CMN"];

nmax=length(DatosCMP)


#=
Ponemos los datos listos para trabajar y les damos un oclayo
=#

datcmp=dictatablaord(DatosCMP)
datcmn=dictatablaord(DatosCMN)

#fig=plot4Ddiscs(datcmp)
#fig=plot4Ddiscs(datcmn)

datcmp=filtraquantil(datcmp)
datcmn=filtraquantil(datcmn)

m4ddistpos=hazmatrizdist(datcmp)
m4ddistneg=hazmatrizdist(datcmn)

herclustpos=hclust(m4ddistpos)
herclustneg=hclust(m4ddistneg)

clustpos=cutree(herclustpos, h=h)
clustneg=cutree(herclustneg, h=h)

tantospos=length(unique(clustpos))
tantosneg=length(unique(clustneg))

println("Obtuvimos $tantospos conjuntos de CM positivos y
        $tantosneg de CM negativos")


#otra visualizacion pa ver como vamos 

scatterclust(datcmp, clustpos)
scatterclust(datcmn, clustneg)


temp=hcat(datcmp, clustpos)
dict4delectrodos=declustaset(temp, enteros=true)
electrodosgruposgrandes=filtraclusterchicos(dict4delectrodos, 3);

for k in keys(electrodosgruposgrandes)
    l=length(electrodosgruposgrandes[k])
    println(" el grupo $k tiene $l electrodos para regiones de Sources")
 end
