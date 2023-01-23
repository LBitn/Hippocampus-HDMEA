module FuncionFranco01

export SigmaVentanas

function SigmaVentanas(Datos, ancho, paso)
    
    ateonde=floor(ancho/2)
    tantos=convert(Int, ceil(length(Datos)/paso))
    result=zeros(tantos)
    for j=ateonde:paso:length(Datos)
        result[j]=std(Datos[j-ateonde:j+ateonde])
        print(j)
        print(result[j])
    end
    return result  
end

end

