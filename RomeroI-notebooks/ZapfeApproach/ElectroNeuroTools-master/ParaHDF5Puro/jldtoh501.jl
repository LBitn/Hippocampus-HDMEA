using HDF5,JLD

archivoname=ARGS[1]
nombregeneral=archivoname[1:end-3]

outname=nombregeneral*"h5"

#las variables que nos sirven

archivo=load(archivoname)

lfp=archivo["LFPSaturados"]
freq=archivo["freq"]*1000 ##en Hz, asi lo quiere Makarov

(alto,ancho,tmax)=size(lfp)

lfp=reshape(lfp,alto*ancho,tmax)


h5open(outname, "w") do file
    write(file, "Args/Fs",freq)
    write(file, "LFP", lfp)
end


