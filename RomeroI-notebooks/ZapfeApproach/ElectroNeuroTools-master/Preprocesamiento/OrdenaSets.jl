module OrdenaSets

push!(LOAD_PATH, "./")
using ArraySetTools
# funciones que ordenan conjuntos respecto a un punto 
# tambien funciones que hacen lineas "rectas" en Z^2


export euclidist, puntoZ2enlinea, barreconjuntoyordena,
sortbydistance, intersectlineamancha, lineaenteros

function euclidist(xxs::Array, yys::Array)
    result=sqrt(sum((xxs.-yys).^2))
end


function puntoZ2enlinea(puntoini::Array, direccion::Number)
    # los puntos estan en "renglon columna"
    # origen abajo izquierda llamado [1,1]
    xo=puntoini[2]
    yo=puntoini[1]
    result=[yo,xo]
    
    if isapprox(direccion, pi, rtol=0.01)
        x1=1
        y1=yo
        result=[y1,x1]
    elseif isapprox(direccion, pi/2, rtol=0.01)
        x1=xo
        y1=64
        result=[y1,x1]
    elseif isapprox(direccion, 0, rtol=0.01)
        x1=64
        y1=yo
        result=[y1,x1]
    elseif  isapprox(direccion, 3*pi/2, rtol=0.01) 
       x1=xo
       y1=1
        result=[y1,x1]
    else
        #cruce pared izquierda:
        xizq=1
        yizq=(xizq-xo)*tan(direccion)+yo
        #cruce pared derecha
        xder=64
        yder=(xder-xo)*tan(direccion)+yo
        #cruce arriba
        yup=64
        xup=xo+(yup-yo)/tan(direccion)
        #cruce abajo
        ydow=1
        xdow=xo+(ydow-yo)/tan(direccion)
        
        if 0<direccion<pi/2
            cruces=[[yder, xder], [yup, xup]]
        elseif pi/2<direccion<pi
            cruces=[[yizq, xizq], [yup, xup]]
        elseif pi<direccion<3*pi/2
            cruces=[[yizq, xizq], [ydow, xdow]]
        else
            cruces=[[yder, xder], [ydow, xdow]]
        end
        
        distan=map(xxs->euclidist(xxs, [yo,xo]), cruces)
        (dd, idx)=findmin(distan)
        
        
        result=cruces[idx]
        result=round.(Int, result)
      #  result=cruces
    end
    
    return result    
        
end

function puntoZ2enlinea(puntoini::Array, puntofin::Array)
    dy=puntofin[1]-puntoini[1]
    dx=puntofin[2]-puntoini[2]
    theta=atand(dy,dx)
    if theta<0
        theta+=2*pi
    end
    result=puntoZ2enlinea(puntoini,theta)
    return result
end
    

function lineaenteros(xx::Array, yy::Array)
    #obtener una lista de pixeles/electrodos a lo largo de una linea con dos endpoints enteros.
    # un punto (renglon columna) y luego el otro
    # asumimos puntos diferentes
    xuno, xdos=xx[2], yy[2]
    yuno, ydos=xx[1],yy[1]
    
    
    
    if abs(xuno-xdos)>= abs(yuno-ydos)
    
        
        longitud=Int(abs(xuno-xdos))+1
        m = (ydos-yuno)/(xdos-xuno)
        
        if xdos>xuno
            rango=xuno:xdos
        else
            rango=xuno:-1:xdos
        end
        xresult=collect(rango) 
        yresult=zeros(Int, longitud)
        for j in 1:longitud
            yresult[j]=round(Int,yuno+(rango[j]-rango[1])*m)
        end
    else
        
        longitud=Int(abs(yuno-ydos))+1
        m=(xdos-xuno)/(ydos-yuno)
        
        if ydos>yuno
            rango=yuno:ydos
        else
            rango=yuno:-1:ydos
        end
        yresult=collect(rango) #asumimos que el ydos es el grande!
        xresult=zeros(Int, longitud)

    for j in 1:longitud
        xresult[j]=round(Int,xuno+(rango[j]-rango[1])*m)
    end
    end
    
    result=hcat(yresult, xresult)
    return result

end
 

function intersectlineamancha(mancha::Array, linea::Array)
    # al parecer no es facil vaciar mancha por la "regla de las referencias "  o algo asi
    lala=rowstoset(linea)
    cha=rowstoset(mancha)
    enlineaymancha=intersect(lala,cha)
   #= if isempty(enlineaymancha)
        print("   eeempt  ")
    end
    =#
    result=elemtorow(enlineaymancha, renglones=true)
   # result=enlineaymancha
    return result
end


function sortbydistance(puntos::Array, referencia::Array)
   fufu(x)=euclidist(x, referencia)
    #permuta los indices para que este el punto mas lejano hasta arriba
   guaga=sortperm(vec(mapslices(fufu, puntos, dims=2)), rev=true)
    # ordena los puntos originales de acuerdo a la permutacion.
    return puntos[guaga, :]
end

function barreconjuntoyordena(datos::Array, punto::Array; dire=pi, dextro=true)

    alfamin=atan(1/64) # minimo paso para cambiar al menos un cuadrito
    if !dextro
        rangodire=dire:alfamin:(dire+2*pi-alfamin) # toda la vuelta levogira
    else
        rangodir=dire:-alfamin:(dire-2*pi+alfamin)
    end

    println("alfamin = ", alfamin)
    result=zeros(Int, 1,2)
    
    for alfa in rangodir
        a=mod2pi(alfa)       
        p=puntoZ2enlinea(punto, a)
        pps=lineaenteros(punto, p)
      #  println(" a = ", a, " p= ", p, "punto= ", punto, " pps= ", pps)
        # print(alfa, " ! ")
        cc=intersectlineamancha(datos, pps)
        
      #  ll=size(cc, 1)
     #   println(" En la direccion alpha=", a, " encontramos ", ll, " puntos.")
        
        if length(cc)>0
            ccoo=sortbydistance(cc, punto)
         
            result=vcat(result, ccoo)
        end
    end
    # no elegante pero funciona (stupid but works)
    result=unique(result[2:end,:], dims=1)
    ll=size(result,1)
    orden=collect(1:ll)
    result=hcat(result, orden)
    return result
    
end

end #module
