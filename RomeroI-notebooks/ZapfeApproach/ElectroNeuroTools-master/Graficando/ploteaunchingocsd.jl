using HDF5
using PyPlot
using Statistics


nombre=ARGS[1]

arxivo=h5open(nombre)
elementos=names(arxivo)
csdlind=read(arxivo["CSDALindenberg"])

sigma=std(csdlind)
mu=mean(csdlind)

function subruido(x,u)
    result=0
    if -u <x < u
        result=0
     else
        result=x
    end
    return result
 end

(alto,ancho,nmax)=size(csdlind)
intervalo=1:5:nmax
#(xmin, xmax)=read(arxivo["cols"])
#(ymin, ymax)=read(arxivo["rengs"])
freq=read(arxivo["freq"])


if length(ARGS) == 1
    nomout="csdplot_"
else
    nomout=ARGS[2]
end


ini=floor(Int, 400*freq)
fin=floor(Int, 500*freq)


(mincsd, maxcsd)=extrema(csdlind)
#limites=min(abs(mincsd),abs(maxcsd))/1.25
limites=5*sigma

close(arxivo)

println(limites, " es la marca de color csd maximo/minimo")

csdadulterada=map(x->subruido(x,2.5*sigma), csdlind) 


csdlind=0
#for interpol in lista


#for n =1:5:nmax
for n=ini:fin

    xlim(0,ancho+1)
    ylim(0, alto+1)
    
    fafa=figure(figsize=(4,4))
    ejemplo=csdadulterada[:,:,n]
    tiempo=round(Int,n/freq)

    title("t = $tiempo ms")
    
   
    imshow(ejemplo, origin="lower",
           cmap="coolwarm", interpolation="bicubic",
           vmin=-limites, vmax=limites,extent=[0.5,ancho+0.5,0.5,alto+0.5])

    nstring=lpad(n, 8, "0")
    
    nomine=nomout*"$nstring.png"

    savefig(nomine, dpi=92)
    close(fafa)
    if mod(n-1,50)==0
        print(" n = ", n, "; ")
    end
end
    
