
#=
Abre y revisa archivos *brw (hdf5) de experimentos con
tejido estriado. Busca canales con posible señal y
canales inservibles segun tres criterios distintos.
Separa el archivo en cachos dependiendo del RAM
disponible. Guarda los resultados como un
hdf5 estandar.
=#

#= Hace escencialmente lo mismo que PreProcEstriado.ipynb pero de
forma no interactiva, automaticamente, pues. Para una versión visual
e interactiva usar el cuaderno de Jupyter =#



push!(LOAD_PATH, ".")

using PreprocTools, AnalisisEstriadoTools

if length(ARGS)==0
    error("Dar el nombre de un archivo *.brw para analizar,
          de preferencia con datos del tejido estriado.")
else
    abrestring=ARGS[1]
    println("Voy a trabajar con el archivo ", abrestring)
end

abrestring=ARGS[1]
stringgeneral=replace(abrestring, ".brw"=>"")  


# comprobamos que tengamos RAM para trabajar
# Arrays muy grandes.
ramlibre=Sys.free_memory()/(1024^3)
if ramlibre>8.1
    tamax=800*1024*1024
else
    tamax=400*1024*1024
end

#establecemos si vamos a partir el archivo en cachos
# y cuantos cachos.
tam=filesize(abrestring)
cachos=div(tam,tamax)+1
println("Voy a dividir tu archivo final en ", cachos, " cachos")

# notas=read(a["3BUserInfo"]["ExpNotes"])
# println(notas) #notas experimentales

Datos=AbreyCheca(abrestring)

freq=Datos["frecuencia"]/1000 #cuadros por milisegundo.
factor=Datos["factor"]

cuadrosmax=Datos["numcuadros"]
tiempototalms=round(cuadrosmax/freq; digits=1) 
println("Tienes ", cuadrosmax, " cuadros de muestreo a ", round(freq; digits=4) ," cuadros por milisegundo")
println( "Esto corresponde a  ", tiempototalms, "ms." )

if size(Datos["DatosCrudos"])[1] != 4096
DatosCrudosArreglados=reshape(Datos["DatosCrudos"], (4096, Datos["numcuadros"]))
else
DatosCrudosArreglados=Datos["DatosCrudos"]
end;


if cachos > 1
    ## escoga el usuario un cacho para trabajar entre 1 y cachos
    cacho=6
    pedazo=round(Int, Datos["numcuadros"]/cachos)
    if cacho < cachos # si no es el ultimo cacho
        desde=1+(cacho-1)*pedazo
        hasta=cacho*pedazo
    else
        desde=1+(cacho-1)*pedazo
        hasta=cuadrosmax   
    end
    
    DatosCrudosArreglados=DatosCrudosArreglados[:,desde:hasta]
    haztodo(DatosCrudosArreglados, freq, factor)
    palabritaespecial="/Cacho_$(cacho)_"
    
else

    haztodo(DatodsCrudosArreglados, freq)
    
    palabritaespecial="/Completo_"

end


println( "Este será el identificador del archivo en cachos: ", palabritaespecial)


#    end;
