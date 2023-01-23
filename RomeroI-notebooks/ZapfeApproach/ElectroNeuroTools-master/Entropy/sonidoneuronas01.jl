using HDF5, StatsBase, Statistics,PyPlot
push!(LOAD_PATH, "../Preprocesamiento/")
using ArraySetTools, OrdenaSets

using WAV

arxname="/home/karel/BRWFiles/Facilitada/Completo_19115s1cut_single_event_preproc.h5"
arx=h5open(arxname)
names(arx)

canales=read(arx["CanalesBuenos"]);
malos=read(arx["CanalesMalos"]);
#lfp=read(arx["dset"]);
lfp=read(arx["LFPSaturados"]);
#freq=read(arx["SamplingRate"])/1000
freq=read(arx["freq"])

dd=size(lfp)
if length(dd)==2
    nmax=dd[2]
    lfp=reshape(lfp, 64,64,nmax);
end
nmax=size(lfp,3)
typeof(lfp)

decibelear=-1
fs=8000

largito=repeat(lfpex, inner=12).*(10.0^decibelear);
wavplay(largito, fs)

lfpexmalo=lfp[62,4,:];

largito2=repeat(lfpexmalo, inner=6).*(10.0^decibelear);
wavplay(largito2, fs)

lfpexrui=lfp[30,30,:];

largito3=repeat(lfpexrui, inner=12).*(10.0^decibelear);
wavplay(largito3, fs)



