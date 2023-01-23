#=
Principal Component Analysis for a Subset of Channels of the BioCAM
=#

using JLD ##Dependiendo que tengamos
using MultivariateStats ##Esto ya trae PCA
using PyPlot ##Just testin'

nombreentrada=ARGS[1]

archivo=load(nombreentrada)

# Variables generales
lfp=archivo["LFPSaturados"]
freq=archivo["freq"]
chans=archivo["Canalesrespuesta"]

(alto,ancho,tmax)=size(lfp)
nchans=length(chans)

lfpselecto=zeros(nchans,tmax)
listacanales=zeros(Int,nchans,3)

for (c, ch) in enumerate(chans)
    listacanales[c,:]=[ch[1],ch[2], ch[1]+64*(ch[2]-1)]
end

listacanales=sortrows(listacanales,by=x->(x[3]))

for c=1:nchans
    lfpselecto[c,:]=lfp[listacanales[c,1],listacanales[c,2],:]
end



##Aplicar PCA, obtener un objeto ideal transformado
## y una reconstruccion del original
tururu=fit(PCA,lfpselecto,maxoutdim=20) 
PCte=transform(tururu,lfpselecto)
rlfp=reconstruct(tururu, PCte)

##Guardemos
archivo["CompPrin"]=PCte
archivo["LFPrecons"]=rlfp
archivo["RespuestaOrden"]=listacanales
save(nombreentrada,archivo)



(npcas,tmax2)=size(Yte)
axis("off")
ylim(-150,1500)
xlim(0,500)
for j=1:npcas
    plot(1:tmax2, Yte[j,:]+(j-1)*100,lw=0.5)
end


axis("off")
for i=1:nchans
    plot(1:tmax,lfpselecto[i,:]+i*100,color="black",lw=0.5)
    plot(1:tmax,rlfp[i,:]+i*100+20, color="red",lw=0.5)
end


yy=lfpselecto-rlfp

for i=1:nchans
    plot(1:tmax,yy[i,:]+i*100+10, color="purple")
end


