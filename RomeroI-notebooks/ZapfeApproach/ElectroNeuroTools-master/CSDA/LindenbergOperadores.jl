module LindenbergOperadores

export GaussianSmooth, DiscreteLaplacian, GaussSuavizarTemporal

#=
Los operadores, funciones y objetos auxiliares que se usan para el CSDA
diferencial (dCSD).
Las funciones auxiliares que lidian sobre conjuntos fueron pasadas
a ParaSets.jl

=#


function UnNormGauss(x,sigma)
    return exp(-x*x/(2*sigma))
end

function GaussSuavizarTemporal(Datos,Sigma=3)
    #Un suavizado Gaussiano temporal.
    #Esto es escencialmente un filtro pasabajos.
    #Depende implicitamente de la frecuencia de muestreo.
    #sigma esta medido en pixeles, es la desviacion estandar de nuestro kernel.
    #El medioancho de nuestra ventana seran 3*sigma

    medioancho=ceil(Sigma*3)
    colchon=ones(medioancho)
    result=zeros(size(Datos))
    datoscolchon=vcat(colchon*Datos[1], Datos, colchon*Datos[end])
    kernel=map(x->UnNormGauss(x,Sigma), collect(-medioancho:medioancho))
    kernel=kernel/(sum(kernel))

    #La convolucion asi normalizada preserva el valor RELATIVO entre los puntos de la funcion.
   
    for t=medioancho+1:length(Datos)+medioancho
        result[t-medioancho]=sum(datoscolchon[t-medioancho:t+medioancho].*kernel)
    end
 
    return result
end


#De momento todo "in file". El Kernel Gaussiano Bidimensional (sigma = 3 pixeles)
GaussianKernel=[0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067
0.00002292	0.00078634	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
0.00038771	0.01330373	0.11098164	0.22508352	0.11098164	0.01330373	0.00038771
0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
0.00002292	0.00078633	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
    0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067]

function GaussianSmooth(Datos)
    tamanodatos=size(Datos)
    result=zeros(tamanodatos)
    temp=copy(Datos)
    (mu, lu)=size(Datos)
    #Primero, hacemos el padding con copia de los datos para que no se suavice demasiado
    ## Okey, parece que los imbeciles de rioarriba cambiaron la sintaxis de
    # rebanadas de matriz. Ahora CUALQUIER rebanada de matriz es colvec.
    arriba=reshape(temp[1,:],(1,lu))
    abajo=reshape(temp[end,:],(1,lu))
    arr3=vcat(arriba,arriba,arriba)
    aba3=vcat(abajo,abajo,abajo)
    
    temp=vcat(arr3, temp, aba3)
    
    for j=1:3
        temp=hcat(temp[:,1], temp, temp[:,end])
    end
    
    for j=4:tamanodatos[1]+3, k=4:tamanodatos[2]+3
        #los indices van primero, "renglones", luego "columnas", etc
        aux=temp[j-3:j+3,k-3:k+3]
        result[j-3,k-3]=sum(GaussianKernel.*aux)
    end
    #Esta convolución no respeta norma L2
    #result=result*maximum(abs(Datos))/maximum(abs(result))
    return result
end


#El operador de Laplace-Lindenberg
LaplacianTerm1=[[0 1 0]; [1 -4 1]; [0 1 0]]
LaplacianTerm2=[[0.5 0 0.5]; [0 -2 0]; [0.5 0 0.5]]
LaplacianKernel=(1-1/3)*LaplacianTerm1+(1/3)*LaplacianTerm2

function DiscreteLaplacian(Datos)
    
    temp=copy(Datos)
    (mu,lu)=size(Datos)
    
    izq=reshape(temp[1,:],(1,lu))
    der=reshape(temp[end,:],(1,lu))
    #Primero, hacemos el padding con copia de los datos para que no se suavice demasiado
    temp=vcat(izq, temp, der)
    temp=hcat(temp[:,1], temp, temp[:,end])
    largo,ancho=size(temp)
    aux=Array{Float32}(undef, 3,3)
    result=zeros(size(temp))
    for j=2:largo-1, k=2:ancho-1
        #los indices van primero, "renglones", luego "columnas", etc
        aux=temp[j-1:j+1,k-1:k+1]
        result[j,k]=sum(LaplacianKernel.*aux)
    end
    #DO  Crop the borders
    result=result[2:end-1,2:end-1]
    return result
end


end


