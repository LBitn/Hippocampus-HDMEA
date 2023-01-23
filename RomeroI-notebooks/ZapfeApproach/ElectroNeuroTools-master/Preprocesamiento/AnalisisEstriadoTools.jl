module AnalisisEstriadoTools

using PreprocTools
using PyPlot, Statistics, HDF5

export haztodo

function haztodo(datacruda::Array, freq::Real, factor::Real,
                 inibu=0.1,
                 finbu=100)
    ## funcion que escencialmente encapsula
    ## las funciones de busqueda de saturados y respuesta.
    ## hace graficas y guarda todo.

    DatosCentrados=FormaMatrizDatosCentrados(datacruda, factor);

    #######################ANALISIS###############################
    
    #busqueda de posibles activos
    PruebaRespuesta=BuscaCanalRespActPot(DatosCentrados,freq,
                                         inibu,finbu,
                                         -85, -300, 11, 24 )
#
# y luego los saturados
    Saturados=BuscaSaturados(DatosCentrados,freq, 800,
                             inibu,finbu)
    
    OtrosSaturados=BuscaSaturadosStd(DatosCentrados, inibu,finbu, freq, 26)

    yOtrosMas=BuscaRuidosos(DatosCentrados, inibu, finbu, freq, 100, 3)
    
    ns=size(Saturados)[1]
    ns2=size(OtrosSaturados)[1]
    ns3=size(yOtrosMas)[1]

    nr=size(PruebaRespuesta)[1]
    println("")
    println("Encontramos ", ns, " canales probablemente saturados.")
    println("Encontramos ", ns2, " canales probablemente con burbuja saturados.")
    println("Encontramos ", ns3, " canales probablemente con ruido positivo.")
    println("Encontramos ",  nr, " canales probablemente con actividad.")


    #######################Inspeccion Visual###############################

    ## Primera Figura: Desviacion Estandar y Canales Marcados
    DesviacionPorCanal=zeros(64,64)
    PromPorCanal=zeros(64,64)
    cini=round(Int, inibu*freq)
    cfin=round(Int, finbu*freq)
    for j=1:64
        for k=1:64
            ChorizoExemplo=vec(DatosCentrados[j,k,cini:cfin])
            DesviacionPorCanal[j,k]=std(ChorizoExemplo)
            PromPorCanal[j,k]=mean(ChorizoExemplo)
        end
    end


    figure(figsize=(7,7))
    xlim(0,65.5)
    ylim(65.5,0)
    title("σ")
    limites=40
    imagen=imshow(DesviacionPorCanal, origin="lower",
                  interpolation="nearest",cmap="inferno", 
                      vmin=0,vmax=limites, extent=[0.5,64.5,0.5,64.5])
    cb=colorbar(fraction=0.046, pad=0.04)

    x=[]
    y=[]
    for j in 1:size(PruebaRespuesta)[1]
        append!(x,PruebaRespuesta[j,2])
        append!(y,PruebaRespuesta[j,1])
    end


    x2=[]
    y2=[]
    for j in 1:size(Saturados)[1]
        append!(x2,Saturados[j,2])
        append!(y2,Saturados[j,1])
    end


    x3=[]
    y3=[]
    for j in 1:size(OtrosSaturados)[1]
        append!(x3,OtrosSaturados[j,2])
        append!(y3,OtrosSaturados[j,1])
    end

    x4=[]
    y4=[]
    for j in 1:size(yOtrosMas)[1]
        append!(x4, yOtrosMas[j,2])
        append!(y4, yOtrosMas[j,1])
    end

    scatter(x,y, marker="x",c="lightblue", s=15)
    scatter(x2,y2, marker="o",c="green", s=15)
    scatter(x3,y3, marker="+",c="green", s=15)
    scatter(x4,y4, marker="|",c="green", s=15)


    savefig("Desviacion.png", dpi=92)


    #### trazos ejemplo

    figure(figsize=(12,20))
    xlabel("ms")
    ylabel("mV")
    ylim(-10,4400)

    d1=34000
    h1=36000
    #intervalo=1:cacho
    intervalo=d1:h1
    intert=intervalo./(freq)


    for j in 1:29
        l=PruebaRespuesta[j,1]
        k=PruebaRespuesta[j,2]
        egtrazo=DatosCentrados[l,k,intervalo].+(j-0.5)*150
        egtrazosuave=mediamov(egtrazo, freq, 0.5)
        annotate("$l, $k ", xy=(1900,(j-0.25)*150), fontsize=14)
        plot(intert,egtrazo, lw=0.5, c="#100070")
        plot(intert,egtrazosuave, lw=3, c="#107710")
    end

    bu=[1776,1796]
    ba=[0,0]

    bru=[1796,1796]
    bra=[0,200]

    #plot(bu,ba, lw=3, c="black")
    #plot(bru,bra, lw=3, c="black")


    #annotate("20ms", xy=(1780,20), fontsize=14)
    #annotate("200mV", xy=(1797,140), fontsize=14, rotation=90)

    #axis("off")
    grid()

    tight_layout()
    savefig("TrazosEjemploEstriado04.png", dpi=92)
    
    return 0
end




end #module

