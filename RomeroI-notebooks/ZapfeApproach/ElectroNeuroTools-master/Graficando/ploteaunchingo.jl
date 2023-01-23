

using JLD,PyPlot

nombre="../JLDFiles/140116s2cut2_evento_4.jld"
#nombre=ARGS[1]

arx=load(nombre)
lfp=arx["LFPTotal"];
freq=arx["freq"];
#csd=arx["kCSDCorrecta"]
(alto,ancho,tmax)=size(lfp)

maxabs=zeros(alto,ancho)

for j=1:64,k=1:64
   maxabs[j,k]=maximum(abs(lfp[j,k,:]))
end

maxmax=maximum(maxabs)

maxabs/=maxmax

#de donde a donde (tiempo) vas a graficar
t1=900
t2=1500
xxs1=[1,t2-t1]
yys1=[400,400]
xxs2=[2,2]
yys2=[-100,400]
intervalo=round(Int,(t2-t1)/freq);

ion()
imagen=figure(figsize=(15,16))
pasox=170
grl=0.075 # ancho de linea
axis("off")
xlim(-10,(t2-t1+paso)*16.1)
ylim(9000,-620)
plot(xxs1,yys1, lw=2,c="black")
plot(xxs2,yys2, lw=2,c="black")
text(50,270,"$intervalo ms")
text(50,20, "0.5mV")
for j=1:4:64 ,k=1:4:64
    #= el evento 8 tiene mal este electrodo
    if (k==57 && j==21 )
        plot((t1:t2)+(k-1)*pasox-t1,-vec(lfp[j+1,k,t1:t2])+(j-1)*125,
         c="black",lw=0.1)
    else
    =#
   # if (k==1) && (j==1) #el experimento 180116 esta al revez
    plot((t1:t2)+(k-1)*pasox-t1,-vec(lfp[j,k,t1:t2])+(j-1)*125,
         c="black",lw=grl)
        #elseif (k==1) && (j==41) #evento 3 del 180116 tiene mal este electrodo
    #    plot((t1:t2)+(k-1)*pasox-t1,-vec(lfp[j+1,64-k,t1:t2])+(j-1)*125,
    #         c="black",lw=grl)
    #else
    #       plot((t1:t2)+(k-1)*pasox-t1,-vec(lfp[j,64-k,t1:t2])+(j-1)*125,
    #     c="black",lw=grl)
    #end
    #end       
end
for j=1:4:64
    text((j+1)*pasox,-300, "$j")
    text(-250,(j-1)*125, "$j")
end
#savefig("MapaTrazosLFPdeshinibida02test.svg",dpi=90)
#savefig("MapaTrazosCSDDeshinibida01.png",dpi=90)
#close()
#close()
