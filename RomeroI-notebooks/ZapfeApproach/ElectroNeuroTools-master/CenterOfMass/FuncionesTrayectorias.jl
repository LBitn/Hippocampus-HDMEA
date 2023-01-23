

function dist2D(x,y)
    #Distancia Euclediana
    result=sqrt((x[1]-y[1])^2+(x[2]-y[2])^2)
    return result
end


function leeunmegaarrayarchivos(NomineGeneralis::AbstractString, desde=35, hasta=501)
    #Cargar una lista de archivos de CM (obsoleto)
    LeMegaArray=Array[]
    for t=desde:hasta
      #  println("$NomineGeneralis-$t.dat")
        CMx=try 
            readdlm("$NomineGeneralis-$t.dat")
        catch
            []
        end
        push!(LeMegaArray, CMx)
    end
    return LeMegaArray
    end


function leunjlddeCM(datos, desde=1, hasta=300)
    # Cargar una lista de CM de un JLD
    LeMegaArray=Array[]
    for t=desde:hasta
      #  println("$NomineGeneralis-$t.dat")
        CMx=try 
        readdlm("$NomineGeneralis-$t.dat")
        catch
            []
        end
        push!(LeMegaArray, CMx)
    end
    return LeMegaArray
end



function encuentraTrayectorias(Datos, mincadena=20)
# Encuentra las trayectorias putativas siguiendo la distancia euclediana

    
    toleradist=8.0*sqrt(2)
    tau=1
    t=1
    j=1
    Catenario=Set{Array{Any}}()
    Cadena=[0 0 0 0]
    CopiaMegaArray=deepcopy(Datos);
    NumFrames=length(Datos)
    FakeNumFrames=NumFrames



    while t <= FakeNumFrames-1 
    
        tau=t
  
        @label arrrrh
 
        if(CopiaMegaArray[tau]==[])     
            jmax,nada=0,0
        else        
            jmax,nada= size(CopiaMegaArray[tau])
        end
    
        
    while j <=jmax && tau<FakeNumFrames
            
            if abs(CopiaMegaArray[tau][j,3]) > 0.05
                
                Eslabon=[transpose(CopiaMegaArray[tau][j,:]) tau]
                Cadena=vcat(Cadena, Eslabon)
         
                mindist=2
                kasterisco=1
                
                if CopiaMegaArray[tau+1]==[]
                    kmax,nada=0,0
                else
                    kmax, nada= size(CopiaMegaArray[tau+1])
                end
                huboalgo=false
          
        
                
                for k=1:kmax
                    
                    EslabonTentativo=CopiaMegaArray[tau+1][k,:]
                    
                    if abs(EslabonTentativo[3])>0.05
                        dist=dist2D(Eslabon,EslabonTentativo)                  
                    if dist<mindist
                        mindist=dist
                        kasterisco=k
                        
                        # println(kasterisco, "=k*", k, "=k")
                        huboalgo=true
                    end
                    end
                
                end    
            
                if huboalgo && mindist<toleradist
                #quitamos el anterior
                    CopiaMegaArray[tau][j,3]=0.0000 
                    # println(mindist," ", t, " ", tau+1 ," ", kasterisco )
                    
                    if tau+1<FakeNumFrames
                        tau+=1
                        j=kasterisco
                        #              println("Pepe t: ", t, "  tau: ", tau, " y  j: ",j )
                        @goto arrrrh
                        
                    else
                        
                        Eslabon=[transpose(CopiaMegaArray[tau+1][kasterisco,:]) tau+1]
                        Cadena=vcat(Cadena, Eslabon)
                        #              println("Pipi t: ", t, "  t: ", t, " y  j: ",j )
                        j+=1
                        tau=t
                        
                        if size(Cadena)[1]>mincadena
                            push!(Catenario, Cadena[2:end,:])
                        end
                        
                        Cadena=[0 0 0 0]
                        @goto arrrrh
                    end
                    
                else
                    
                    if size(Cadena)[1]>mincadena
                        push!(Catenario, Cadena[2:end,:])
                    end
                    Cadena=[0 0 0 0]
                    j+=1
                    tau=t
                    @goto arrrrh
                end
                
            end #cierra sobre el if de  la masa 
        
        j+=1                    
        tau=t
        
    end
    @label urrr
            
     j=1
     t+=1
     tau=t
     Cadena=[0 0 0 0]
     end 
        
    return Catenario

 end
    


