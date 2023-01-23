using JLD , PyPlot


if length(ARGS)==0
    error("Dar el nombre de un archivo LFP_*.jld para graficar")
else
    archivo=ARGS[1]
    println("Voy a trabajar con el archivo ", archivo)
end
    

#Los datos para fraficar
LFP=load(archivo)["LFPSaturados"]
saturados=load(archivo)["CanalesSaturados"]
respuestas=load(archivo)["Canalesrespuesta"];
frecuencia=load(archivo)["freq"]
retrazo=load(archivo)["retrazo"]


palabra=replace(archivo,".jld", "")
println(palabra)
palabrita=replace(palabra,"LFP_Promedio_", "")

if !(isdir(palabrita))
     mkdir(palabrita)
end

#parametros
tantossaturados=length(saturados)
tantosrespuesta=length(respuestas);
limites=200

xxsresp=zeros(tantosrespuesta)
yysresp=zeros(tantosrespuesta)
j=1
for q in respuestas
    xxsresp[j]=q[2]
    yysresp[j]=q[1]
    j+=1
end
xxssatu=zeros(tantossaturados)
yyssatu=zeros(tantossaturados)
j=1
for q in saturados
    xxssatu[j]=q[2]
    yyssatu[j]=q[1]
    j+=1
end


ioff()
tmax=size(LFP,3)
#tmax=4 #for tests

#= El mapa de color del potencial y los canales
saturados y posibles respuestas
=#




for t=1:tmax
   
    figure(figsize=(6,6))
    xlim(0,65)
    ylim(0,65)

    ejemplo=LFP[:,:,t]
    tiempo=round((t-retrazo)/frecuencia,1)
    title("t= $tiempo")
    imagen=imshow(ejemplo, origin="lower", interpolation="nearest", cmap="Spectral_r", 
vmin=-limites,vmax=limites, extent=[0.5,64.5,0.5,64.5])
    cb=colorbar(imagen, fraction=0.046)
    if t>(tmax-140)

        scatter(xxsresp,yysresp, alpha=0.9, c="red", edgecolor="none")
        scatter(xxssatu,yyssatu, alpha=1.0, c="white", edgecolor="white", marker="x", lw=1) 
    end
    nombre=palabrita*"/"*palabra*"_"*string(t)*".png"
    savefig(nombre,dpi=90)
    close("all")
  #  PyPlot.close_queued_figs()
end



#=
Los trazos de las posibles respuestas
=#
desde=round(Int, ceil(5*frecuencia))
tiemporeal=map(x->(x-desde)/frecuencia,1:tmax)
limites=4000

figure(figsize=(8,6))
xlim(-6,12)

ylim(-limites,limites)
thewholefockenevent=[]
promedio=zeros(tmax)
k=0

for q in respuestas
    x=q[2]
    y=q[1]
    if x != 1 && y != 1 
    k+=1
    datos=vec(LFP[y,x,1:tmax])
        thewholefockenevent=vcat(thewholefockenevent,datos)
    promedio+=datos
 
        plot(tiemporeal, datos, c="grey", alpha=0.5)   
    end
end
promedio=promedio/k
plot(tiemporeal, promedio, c="black")
xlabel("tiempo [ms]")
ylabel("LFP [µV]")

tight_layout()
println(k)
outname=string(palabrita,"/","Trazos_",palabra,"-01.png")

savefig(outname, dpi=90)
close()


#=

Los trazos de los canales saturados
=#

figure(figsize=(8,6))
xlim(-6,12)
ylim(-limites,limites)
promedio=zeros(tmax)
k=0
anotherwholefockenevent=[]
for q in saturados
    x=q[2]
    y=q[1]
    if x != 1 && y != 1 && x>45 
    k+=1
    datos=vec(LFP[y,x,1:tmax])
        anotherwholefockenevent=vcat(anotherwholefockenevent,datos)
    promedio+=datos
 
        plot(tiemporeal, datos, c="grey", alpha=0.5)   
    end
end
promedio=promedio/k
plot(tiemporeal, promedio, c="black")
xlabel("tiempo [ms]")
ylabel("LFP [µV]")

tight_layout()
println(k)
outname=string(palabrita,"/","Saturados_",palabra,"-01.png")

savefig(outname, dpi=90)
close()


