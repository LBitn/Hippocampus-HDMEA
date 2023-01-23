function OtsuMethod(Data)
    #Pa que esto funcione hay que mandar los datos ya sin NaNs u otras mugres   
    binsdefault=2*ceil(Int,sqrt(length(Data)))
    (rango, cuentas)=hist(Data,binsdefault)
    tantos=length(rango)
    valores
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

