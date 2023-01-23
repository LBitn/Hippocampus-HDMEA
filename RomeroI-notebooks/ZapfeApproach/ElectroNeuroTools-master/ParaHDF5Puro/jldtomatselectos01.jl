using MAT,JLD


archivname=ARGS[1]
archivo=load(archivname)

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


outname=archivname[1:end-3]*"mat"

outfile=matopen(outname,"w")
write(outfile, "lfp", lfpselecto)
write(outfile, "channels", listacanales)
write(outfile, "Fs", freq*1000)
close(outfile)



texg=137


# Hagamos un mapa de electrodos selectos
estandar=std(lfp[:,:,texg:end],3)
estandar=reshape(estandar, (64,64))

using PyPlot
figure(figsize=(8,8))
xlim(0,65)
ylim(0,65)
grid()
imagen=imshow(estandar, origin="lower", interpolation="nearest", cmap="viridis", 
vmin=0,vmax=24, extent=[0.5,64.5,0.5,64.5])
cb=colorbar(imagen, fraction=0.046)
scatter(listacanales[:,2],listacanales[:,1],s=7, alpha=0.9, c="red", edgecolor="none", label="Respuestas")   
legend()
savefig("MapaRespuestas_v03.png",dpi=90)
close()


