##Cambiamos las convenciones a Julia 0.5 
##En particular, hist e hist! se cambian por fit(Histogram...) en StatsBase
### espero que no conflictue despues con MultivariateStats (no creo).

using StatsBase 

function OtsuMethod(Data)
    #Pa que esto funcione hay que mandar los datos ya sin NaNs u otras mugres   
    binsdefault=2*ceil(Int,sqrt(length(Data)))
    histograma=fit(Histogram,Data,nbins=binsdefault)
    ##por alguna razon los devs creen que esta padre poner el rango en una tupla
    (rango, cuentas)=(histograma.edges[1],histograma.weights)
    tantos=length(rango)
    #valores
    omega1=0
    omega2=0
    mu1=0
    mu2=0
    sigmab=0
    sigmabtemp=0
    tbest=0
    varlim=0
    for t=1:tantos-1
        omega1=sum(cuentas[1:t])
        omega2=sum(cuentas[t+1:tantos-1])
        mu1=sum(cuentas[1:t].*rango[1:t])/omega1
        mu2=sum(cuentas[t+1:tantos-1].*rango[t+1:tantos-1])/omega2        
        sigmabtemp=omega1*omega2*((mu1-mu2)^2)
        if sigmabtemp>sigmab
            sigmab=sigmabtemp
            tbest=t
            varlim=rango[t]
        end
    end
    return (sigmab,tbest,varlim)
end

function OtsuUmbralizar(DatosMatriz)
    # aplanar datos
    DataFlatten=reshape(DatosMatriz, size(DatosMatriz)[1]*size(DatosMatriz)[2])
    mascara=zeros(DatosMatriz)
    umbral=OtsuMethod(DataFlatten)[3]
    mascara=map(x->(x>umbral)?1:0, DatosMatriz)
    return mascara
end

function ParchaMascara01(Datos)
    result=zeros(Datos)
    temp=copy(Datos)
    (alto,ancho)=size(Datos)
    #Primero, hacemos el padding de los datos para que no se suavice demasiado
    
    colchonvertical=zeros(1,ancho)
    colchonhorizontal=zeros(alto+2)
    
    temp=vcat(colchonvertical, temp, colchonvertical)
    temp=hcat(colchonhorizontal, temp, colchonhorizontal)
    
    
    for j=2:alto+1, k=2:ancho+1
        aux=temp[j-1:j+1,k-1:k+1]
        test=sum(aux)
        if(temp[j,k]<2.0 && test<2.1)
            result[j-1,k-1]=0
        #    println(j,k)
        elseif(temp[j,k]<1.1 && test>2.0)
            result[j-1,k-1]=1
        else
            result[j-1,k-1]=temp[j,k]
        end
    end   

    return result
    
end
