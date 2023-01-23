#=
Dividamos un brw enorme en cachos manejables en memoria
Basado en el corta datos de Isabel.
=#

push!( LOAD_PATH, "." ); 
using Z_auxiliaresBRW, Z_auxiliaresGEN
using HDF5

arxname=ARGS[1]

#memoria usada
ramsiz=Sys.total_memory()
pmem = filesize(arxname)/ramsiz


#000000000000 Funciones chidas

println( "Hola Carnal. Si tu archivo mide mas de un octavo del ram
         disponible lo voy a partir en 8 cachos o mas. Si no no hare nada.")

println(" Tu memoria usada al cargar el archivo resulta en $pmem del total.")


function guardadictydatos(oldie::Dict, nuevos::Array, nombre::String)
    #=
    Funcion que guarda un diccionario obtenido de un brw, pero
    remplazando los datos crudos por el array "nuevos"
    =#
    
    llaves=keys(oldie) 
    if in( "dset", llaves)
        
        oldie["dset"]=nuevos
        println("guardando en $nombre")
        
        h5open(nombre, "w")  do file
            for k in llaves
                datos=oldie[k]
                write(file, "$k", datos)
            end
        end
        return
    else   
        error("Tu dict no trae dataset a remplazar")
    end
        
end


function cortayguarda(filenombre::String, pmem=0.125)
    
    
    if pmem<0.125 #pmem es memoria usada, no libre!!
        println("tu archivo tiene tamaño razonable, dejalo asi")
        return 0
    else

    dicth5 = brw_things( filenombre ); # variables útiles
    data = dicth5[ "dset" ]; # path al dataset (NO READ! to much)
    M=size(data, 1) # el tamaño de data
    chans=dicth5["Chs"]
    nch=size(chans, 1)
    nframes=Int(M/nch)
    

    nrams=ceil(Int, filesize(filenombre)*8/Sys.total_memory())  
    ncachos=max(nrams, 8) # al menos partela en 4
    fcacho=floor(Int, nframes/ncachos) # frames por cacho dara un poco menos que el total.
        
    paso=nch*fcacho #cuantos datos tiene cada cacho.
        
    println("Cuadros totales: \t", nframes, "\n", 
            "Datos totales: \t",  M, "\n", 
            "Cachos: \t " ,ncachos, "\n",
            "Cuadros por cacho: \t", fcacho, "\n")
    
    
    for n = 1:ncachos  # numero de Β a cortar (1-> 4096,1->ω)
        
        cini=(n-1)*fcacho+1
        cfin=n*fcacho
        datitos =reshape(data[((n-1)*paso+1):n*paso], nch, fcacho)
        
        nstring=lpad(n, 2, "0")
        dir=dirname(filenombre)
        fn=basename(filenombre)
        nuname=dir*"/Cacho_"*nstring*"_"*fn
        
            
            
        notacacho="Este es un intervalo de un dataset mas grande. \n
                    Contiene los cuadros de $cini a $cfin."
        dicth5["notacacho"]=notacacho
            
        guardadictydatos(dicth5,datitos, nuname)
        
        println( "vamos en el cacho", n,  " de ", ncachos, " que va del cuadro ",
                cini, " al ",  n*cfin, " guardado en $nuname")
    end
    print("bla")
    return  "Chido"
    end

end

#000000000000 End Funciones chidas

# Lo merito bueno

cortayguarda(arxname, pmem)

