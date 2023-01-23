using JLD
include("FuncionesTrayectorias.jl")


numarg=length(ARGS)
if numarg<1
    error("Dame un nombre de archivo JLD que tenga CM")
else
    nombre=ARGS[1] 
end

nombregeneral=nombre[1:end-4]

nombresalida=string(nombregeneral, "-Tray.jld")

println("este es el nombre salida: ", nombresalida)

if numarg==1
    inicadena=0
    fincadena=400
elseif numarg==2
    inicadena=0
    fincadena=ARGS[2]
elseif numarg==3
    inicadena=ARGS[2]
    fincadena=ARGS[3]
end

    


archivo=load(nombre)
DatosCMP=archivo["CMP"]
DatosCMN=archivo["CMN"]

tolerancia=3

CatenarioPositivo=encuentraTrayectorias(DatosCMP,tolerancia,inicadena,fincadena);
CatenarioNegativo=encuentraTrayectorias(DatosCMN,tolerancia,inicadena,fincadena);

println("procedo a escribir en disco duro")
save(nombresalida, "CatenarioNegativo", CatenarioNegativo, "CatenarioPositivo", CatenarioPositivo)
println("cha chan")
