## Preprocesamiento de un LFP para implementar ICA
using JLD

##Funciones de Otsu
include("./Otsu05.jl")

#= Esto es muy importante:
===========================
Si usas OtsuUmbral sobre un LFP con canales Saturados
adivina que pasa?
Exacto: te da los canales Saturados.
=#



archivo=ARGS[1]
println("cargando lfp")
lfp=load(archivo)["LFPSaturados"]

(ancho,alto,tmax)=size(lfp)
### para que sea mas facilito aplastamos los datos de entrada en una lista larga

#probando
(ancho,alto,tmax)=(64,64,2000)

lfp=randn(ancho,alto,tmax)

aa=ancho*alto


lfpplano=reshape(lfp, (aa,tmax))


### Primero los centramos

lfpmasajeado=zeros(lfpplano)
promedios=zeros(aa)
mean!(promedios,lfpplano)
for j=1:aa
    lfpmasajeado[j,:]=lfpplano[j,:]-promedios[j]
end

correlaciones=zeros(aa,aa)

vcdot(x,y) = dot(vec(x), vec(y))

vcdot(lfp[1,:],lfp[1,:])

##Al parecer no podemos usar cov de julia porque normaliza
@time for j=21:aa,k=1:aa
    correlaciones[k,j]=vcdot(lfpplano[k,:], lfpplano[j,:])
end
    

correlaciones/=tmax

writedlm( "test.dat" ,correlaciones)

