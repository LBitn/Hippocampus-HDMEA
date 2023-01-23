module EntropySymb

#=
Modulo para hacer simbólica dinámica de tres simbolos.
para el CSDA de señales electrofisiologicas
=#

using Base.Iterators

export simboliza, problock, ncartprod, issubvec, allprobsn, shanonH_n

epsi=21

function simboliza(x, delta=epsi)
    # funcion que pasa floats a tres simbolos.
    result=0
    if x>delta
        result=1
    elseif x< -delta
        result=-1
    else
        result=0
    end
    return result
end



function problock(xxs::Array, yys::Array)
    # calculamos la prob exp de encontrar el array yys en el array xxs
    lx=length(xxs)
    ly=length(yys)
    lp=lx-ly+1
    result=0
    if(lx<ly)   
     error("La secuencia objetivo es mayor a la secuencia para buscar.")
    else
        for j=1:lp
            if xxs[j:j+ly-1]==yys
                result+=1
            end
        end
    end
    return result/lp
end


function ncartprod(ss, n)
    # el producto cartesiano de los simbolos ss de n dimensiones.        
    argumento=ntuple(i->ss, n)
    result=product(argumento...)
    return(result)
end

# con permiso de StackOverflow oooo.key takes lot to compile!
issubvec(v,big) = any([v == big[i:(i+length(v)-1)] for i=1:(length(big)-length(v)+1)])

function allprobsn(xxs, ss, n)
    # calcula las probabilidades mayores a cero de encontrar cada
    # secuencia de longitud n en ss a partir de los simbolos ss
    # devuelve las no encontradas aparte.
    result=Dict{Array, Float64}()
    palabras=ncartprod(ss, n)
    malas=Set()
    for p in palabras
        pr=[p...]
        aux=problock(xxs, pr)
        if aux>0
            result[pr]= aux
        else
            push!(malas, pr)
        end
    end
    return (result, malas)
end


function allprobsn(xxs, ss, n, vacias)   
    # Esta version puede usar un conjunto de secuencias
    # de longitud m<n que ya no habria que buscar.
    # cuando m es pequeño como 3 o 4 es mas eficiente. 
    result=Dict{Array, Float64}()
    palabras=ncartprod(ss, n)
    malas=Set()
    for p in palabras
        pr=[p...]
        tanteada=any(issubvec(v, pr) for v in vacias)
        if tanteada
         #   print(pr, " esta tanteada!!  ")
            aux=0
        else
            aux=problock(xxs, pr)
            if aux>0
                result[pr]= aux
            else
                push!(malas, pr)
            end
        end # sobre tanteada
    end #sobre palabras
    return (result, malas)
end


    function shanonH_n(xxs, ss, n)
    (probs, v)=allprobsn(xxs,ss,n)
    h=0
    for q in values(probs)
        if q>0.00001
            h+=q*log(q)
        end
    end
    result=-h/n
    return result
end

function shanonH_n(xxs, ss, n, v)
    (probs, v)=allprobsn(xxs,ss,n, v)
    h=0
    for q in values(probs)
        if q>0.00001
            h+=q*log(q)
        end
    end
    result=-h/n
    return result
end

function shanonH_n(probs,n)
    #meter un dict con probs>0 porfis
    h=0
    for q in values(probs)
        h+=q*log(q)
    end
    result=-h/n
    return result
end
    
end #module

