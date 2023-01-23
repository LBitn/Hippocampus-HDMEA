#= Detección de canales activos por 1. metodo de Otsu =#

push!(LOAD_PATH,".") # folder que continen los modulos. Correr desde ahí dentro.
# LOad Modules
using Z_auxiliaresBRW, Z_auxiliaresGEN
# Packages
using JLD, Statistics

#
Γ = OS( "linux" );
# donde estan los jld para procesar
BRWname = "/home/isabel/Desktop/200819/control_01.brw"; # original BRW file
PROCESSpath = string( split( BRWname, "." )[ 1 ], Γ, "Process" );
υ = brw_things( BRWname );
SamplingRate = υ[ "SamplingRate" ];
ξ = υ[ "Chs" ];
# donde se guardarán los archivos de salida de los procesos de detección (de haber)
WORKpath = string( split( PROCESSpath, "Process" )[ 1 ], "Detection", Γ);
checkpath( WORKpath ); # checar si el folder para guardar los bines existe
#
files = find_files( PROCESSpath, ".jld" ); # files ya con path.
#

#= 1. Método de Otsu
1.1 Con ventanas de Desviación Estandar =#
data = load( files[ 1 ] )[ "data" ]; # carga previa, solo para sacar los valores necesarios
# divisores del numero de frames por bin
δφ = div_n_ab( size( data, 2 ), 1, size( data, 2 ) ); 
# mismos divisores pasados a msec (por dimensionar)
δτ = ( δφ./SamplingRate )*1000; 
δφ = δφ[ 1 .<= δτ .<= 500 ]; # seleccionar valores entre 1 y 500 msec

if length( δφ ) <= 10
    ω = δφ; # tamaño final de las ventanas de Otsu a usar
else
    ω = δφ[ iseven.( collect( 1:1:length( δφ ) ) ) ];
end

#= espacio para guardar los valores de STD de cada diferente ventana, 
ej. 5 ventanas, para 65 archivos = 325 columnas, 
con 4095 filas por cada canal menos el de referencia =#
Π = zeros( Int, size( ξ, 1 ) - 1, size(  ω, 1 )*( size( files, 1 ) ) );

for i = 1:size( files, 1 );
    # para cada cacho, se hace una ventana de tiempo previamente definida
    β = zeros( Int, size( ξ, 1 ) - 1, size(  ω, 1 ) );
    #
    data = load( files[ i ])[ "data" ];
    data = data[ setdiff( 1:end, 1 ), : ]; # remueve el canal de referencia (1,1)
    for j = 1:size(  ω, 1 ) # para cada ventana de tiempo
        # limites para tomar cada cacho de la ventana de tiempo
        lims = collect( 1: ω[ j ]:size( data, 2 ) );
        # aqui se colectarán las stds de cada bin del cacho
        α = zeros( size( ξ, 1 ) - 1,size( lims, 1 ) - 1 ); 
        for k = 1:size( lims, 1 ) - 1 # por cada bin del cacho
            # std de cada bin
            α[ :, k ] = std( data[ :, lims[ k ]:lims[ k + 1 ] ], dims = 2 ); 
        end
        α[ α .== Inf ] .= 0; # sometimes the std gives back Inf, du no why
        # aplicando otsu al vector de stds maximas
        β[ :, j ] = OtsuThr( maximum( α, dims = 2 ) ); 
    end
    Π[ :, 1 + ( i - 1 )*size(  ω, 1 ):i*size(  ω, 1 ) ] = β;
    println( "file ", i, " obtained" );
end
save( string( WORKpath, "Otsu_STD_windows.jld" ), "data", Π , "windows", ω );
println( "Otsu method for std windows saved" )


