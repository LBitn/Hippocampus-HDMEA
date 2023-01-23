module Z_auxiliaresGEN
#-------------------------#

function div_n_ab( n::Int, lo::Int, hi::Int )
#=
divisores del numero n entre los rangos lo and hi (Define el numero de cacho en los que quieres cortar el archivo total)
=#
    ρ = collect( 1:Int( floor( sqrt( n ) ) ) ); # los numeros de 1 en 1 de la raiz cuadrada 
    σ1 = findall( n.%ρ .== 0 ); # divisores de la raiz cuadrada ( residuo = 0 )
    σ2 = Int.( ( n )./( σ1 ) ); # Sacar los pares ( de 100, 2-50, 10-10, etc)
    σ = sort( unique( vcat( σ1, σ2 ) ) ); # remover duplicados, concatenar, ordenar
    if @isdefined lo
        rn = σ[ findall( hi .>= σ .>= lo ) ];
        println( "Los divisores de ", n, " entre ", lo, " y ", hi, " son ----- ", rn, " ----- ", size( rn, 1 ), " divisores." );
        return rn
    else
        return σ
        println( "Los divisores de ", n, " son ------ ", σ, " ----- ", size( σ, 1 )," divisores." );
    end
end
#-#
function add_zeros( n::Int, i::Int )
#=
añade el numero n de 0 al valor i de ser necesario
=#
    bin_n = string( i );
    for a = 1:n - 1
        if i < 10^a
            bin_n = string( "0", bin_n );
        end
    end
    return bin_n
end
#-#
function OS( OS::String )
# para los paths de mi lap o del lab

    if OS == "linux"
        com = "/"
    elseif OS == "windows"
        com = "\\"
    end
    return com
end
#-#
function checkpath( workpath::String )
# checar si el folder existe, si no, hacerlo
    if isdir( workpath ) == false
        mkpath( workpath )
    end
end
#-#
function find_files( path_main::String, key::String )
# buscar en el directorio Path_main, los archivos que terminen con key
	searchdir( path_main::String, key::String ) = filter( x -> endswith( x, key ), readdir( path_main ) );
	files = searchdir( path_main, key );
    if path_main[ end ] == '/'
        files = string.( path_main, files );
    else
        files = string.( path_main, "/", files );
    end
	return files
end
#-#
function neighborgs( center::Int, d::Int )
    #=
    obtienen la "d"-vecindad del canal "center"
    =#
    A = reshape( 1:4096, 64, 64 );
    A = Int.( A' );
    x_c = findall( A .== center )[ ][ 2 ]
    y_c = findall( A .== center )[ ][ 1 ]
    aux = [ ( x_c - d ),( x_c + d ), ( y_c - d ), ( y_c + d ) ]
    aux[ aux .< 1 ] .= 1;
    aux[ aux .> 64 ] .= 64;
    neigh = A[ aux[ 3 ]:aux[ 4 ], aux[ 1 ]:aux[ 2 ] ];
    return neigh
end
#-#
using StatsBase
function OtsuMethod(datos::Array)
    #=
    Metodo de Otsu para separar una imagen en dos subconjuntos.
    =#
    min = minimum( datos ); max = maximum( datos );
    #
    binsdefault = 2*ceil( Int, sqrt( length( datos ) ) ); # <-?
    #
    h = fit( Histogram, vec( datos ), closed = :right, nbins = binsdefault );
    ( Y, X ) = ( h.edges[ 1 ], h.weights ); # stupid Statistics, dice Karel
    T = length( Y );
    # valores
    ω1 = 0; ω2 = 0;
    μ1 = 0; μ2 = 0;
    σ = 0; σtemp = 0;
    τbest = 0; varlim = 0;
    #
    for τ = 1:T - 1
        ω1 = sum( X[ 1:τ ] );
        ω2 = sum( X[ τ + 1:T - 1 ]);
        μ1 = sum( X[ 1:τ ].*Y[ 1:τ ] )/ω1;
        μ2 = sum( X[ τ + 1:T - 1 ].*Y[ τ + 1:T - 1 ] )/ω2;
        σtemp = ω1*ω2*( ( μ1 - μ2 )^2 );
        if σtemp > σ
            σ = σtemp;
            τbest = τ;
            varlim = Y[ τ ];
        end
    end
    return ( σ, τbest, varlim )
end
#-#
function OtsuThr( datos::Array ) # depende de OtsuMethod definido anteriormente
    # aplanar datos
    result = zeros( size( datos ) );
    thr = OtsuMethod( datos )[ 3 ];
    result = map( x -> ( x > thr ) ? 1 : 0, datos );
    return result
end

#----#
export div_n_ab, add_zeros, OS, checkpath, find_files, neighborgs, OtsuThr

#-#
end
