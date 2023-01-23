module FFTWtoolkit

using FFTW

export filtbanda, fftfiltrobando, separabandaneuro, BandaNeuro
export PeriodoGramNeuro, periodogramneuro
export MatrizDatoBanda, sacalasbandas
export MiFFTWindow, fftwwin

function filtbanda(datosft::Array, freqs::Any, bajo, alto)
    
    xx=deepcopy(datosft)
    
    ene=length(xx)
    bix=findfirst(x->x>=bajo, freqs)
    aix=findfirst(x->x>alto, freqs)
  
    # buscar el centro de la T. Fourier  
    enemedios=div(ene,2)+1
    
    # matar todo fuera de la banda
    for j=1:bix
        xx[j]=0
        if j>=2
        xx[ene-j+2]=0
        end
    end
    
    for j=aix:enemedios
        xx[j]=0
        xx[ene-j+2]=0
    end
    
    return xx

end
    
function fftfiltrobanda(datostiempo::Array, bajo::Real, alto::Real;  fs=7022)
    #Filtra usando fftw, a lo GUEY.
    if bajo>alto
        alto, bajo = bajo, alto
    end
    
    ene=length(datostiempo)
    fftdatos=fft(datostiempo)
    
    freqs=FFTW.fftfreq(ene, fs)
   
    result= ifft(filtbanda(fftdatos, freqs, bajo, alto))
    
    
end


struct MiFFTWindow
    freqs::Array{Number}
    tiemps::Array{Number}
    valores ::Array{Number}
end

function fftwwin(datos::Array, ventana; fs=7022)
    # la ventana en indices.
    if typeof(ventana)<:Integer
        idn=ventana
    else
        idn=floor(Int64, freq/ventana)
    end
    
    ene=size(datos,1)
    freqs=collect(fftshift(FFTW.fftfreq(idn, fs)))
    plano=plan_fft(datos[1:idn])
    n=ene-idn
    tiempos=collect(((0:(n-1)).+0.5)./fs)
    
    result=zeros(Complex, idn, n)
    
    for j=1:n
        # println(" $j $n $idn")
        result[:,j]=fftshift((plano*datos[j:j+idn-1]))
       # println(" $j $n $idn")
    end
    MiFFTWindow(freqs,tiempos,result)
end



struct BandaNeuro
    signal :: Array
    freqs :: Array{Number}
    names :: Array{String}
end    


function separabandaneuro(datos::Array; fs=7022, planes=(nothing, nothing))
    
    if planes==(nothing, nothing)
        fdatos=fft(datos)
    else
        fdatos=planes[1]*datos
    end

    ene=size(datos,1)
    freqs=FFTW.fftfreq(ene, fs)
    
    result=zeros(ene,7)
    
    kys=["Delta", "Theta", "Alpha", "Beta", "Gamma", "spr", "uspr"]
    frehz=[3.5, 7, 15, 35, 100, 250, 600]
   
    aux=filtbanda(fdatos, freqs, 0, frehz[1])
  
    if planes==(nothing, nothing)
        result[:,1]=real(ifft(aux))
    else
        result[:,1]=real(planes[2]*(aux))
    end
        
    
    for j=2:7
        aux=filtbanda(fdatos, freqs, frehz[j-1], frehz[j])
      
        if planes==(nothing, nothing)
            result[:,j]=real(ifft(aux))
        else
            result[:,j]=real(planes[2]*(aux))
        end
    end
    BandaNeuro(result, frehz, kys)
end
    

struct PeriodoGramNeuro
    powers :: Array{Number}
    freqs :: Array{Number}
    names :: Array{String}
end

   
function periodogramforneuro(datos::Array; fs=7022)
  
    kys=["Delta", "Theta", "Alpha", "Beta", "Gamma", "spr", "uspr"]
    frehz=[3.5, 7, 15, 35, 100, 250, 600]
    pows=zeros(7)
    pp=periodogram(datos, fs=fs)

    idx=findfirst(x->x>frehz[1], pp.freq)
    pows[1]=sum(pp.power[1:idx])
    
    for k=2:7     
        idx1=findfirst( x->frehz[k-1] >= x, pp.freq)
        idx2=findfirst( x-> x >frehz[k], pp.freq)
      
        pows[k]=sum(pp.power[idx1:idx2])
    end
 
    return PeriodoGramNeuro(pows,frehz, kys)
end


struct MatrizDatoBanda
    delta::Array
    theta::Array
    alpha::Array
    beta::Array
    gamma::Array
    spr::Array
    uspr::Array
end

function sacalasbandas(datos::Array; fs=7022)
    tamano=size(datos)
    resultillo=zeros(tamano...,7)

    planida=plan_fft(datos[4,4,:])
    planvuelta=plan_ifft(datos[4,4,:])
    
    for j=1:64, k=1:64
        aux=separabandaneuro(datos[j,k,:],fs=fs, planes=(planida,planvuelta))
        for l=1:7
            # la l-esima banda va en la 4to indice igual a l 
            resultillo[j,k,:,l]=aux.signal[:,l]
        end
    end
    
    MatrizDatoBanda(resultillo[:,:,:,1],resultillo[:,:,:,2],
        resultillo[:,:,:,3],resultillo[:,:,:,4], 
    resultillo[:,:,:,5],resultillo[:,:,:,6],resultillo[:,:,:,7])
    
end


### Probablemente Obsoloeto
function wfft(datos::Array, cuadros::Int)
#hace un fftw por ventanas, lo que hace que tengamos un array de una dimension mas.
    r=1:cuadros
    nredux=length(datos)-cuadros
    plan=plan_fft(datos[r])
    result=zeros(Complex, nredux, cuadros)
    for j=1:nredux
        raux=r.+j
        result[j,:]=plan*datos[raux]
    end
    return result
end


end #module


