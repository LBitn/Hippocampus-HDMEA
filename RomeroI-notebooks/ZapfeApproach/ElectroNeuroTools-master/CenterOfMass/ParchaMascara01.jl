putama=readdlm("MascaraOtsu01.dat")

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

kuak=ParchaMascara01(putama)
koraka=putama-kuak

    
writedlm("Parchada01.dat", kuak)
writedlm("Diff.dat", koraka)


cuaracua=ParchaMascara01(kuak)

writedlm("Parchada02.dat", cuaracua)
    #I like this

    
    
