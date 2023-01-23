#Funciones que encuentran los golpes evocadas y encuentran 

function EncuentraTrancazosRaw(datos::Array, tolerancia=1400)
    #=Esta funcion esta hecha pensando en que se le alimenta
    Raw data, es decir, los enteros que escupe el BioCAM.
    Para usar datos en microvolts es necesario otra funcion
    La función toma el arreglo y busca los datos que tengan
    deflecciones del promedio mayores a tolerancia.
    =#
    result=Int[]
    longitud=length(datos)
    jcomp=0
    media=mean(datos)
    for j=1:longitud
        if abs(datos[j]-media)>tolerancia
            if j-jcomp>1
                push!(result,j)
            end
            jcomp=j
        end
    end
    return result
end




function ActivAlrededorTrancazo(Lista::Array, xxs::Array)
    #= Funcion que toma un pequeño subarray alrededor de los
puntos temporales que selecciona EncuentraTrancazosRaw (u otra similar),
llamados Lista. De ahí se toman 36 cuadros hacia atras (aprox 5ms) y 
211 hacia adelante (aprox 30 ms).  y se devuelve un diccionario con los arrays
=#

    #Aqui no se le ha hecho reshape a las matrices todavia
    result=Dict{AbstractString, Array}()
    q=1
    desde=36
    hasta=211
    for j in Lista
        nomineclave="Trancazo_$q"
        result[nomineclave]=xxs[:,j-desde:j+hasta]
        #println(nomineclave)
        q+=1
    end
    return result
end

function ActividadFueraTrancazo(Lista::Array, xxs::Array)
    #=Lo opuesto de la anterior. Devuelve un array con todo
    menos lo que se encuentre alrededor de los trancazos.
    Buena para hacer estadísica sobre el ruido y la actividad
    de fondo. =#
    
    q=1
    desde=36
    hasta=211
    aux=trues(xxs)
    for j in Lista
        aux[j-desde:j+hasta]=false
    end
    result=zeros(1)
    aux2=find(aux)
    for j in aux2
        result=vcat(result,xxs[j])
    end
    return result
end
        
    
