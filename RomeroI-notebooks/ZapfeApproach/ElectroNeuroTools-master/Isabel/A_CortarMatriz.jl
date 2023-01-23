push!( LOAD_PATH, "." ); # folder que continen los modulos. Correr desde ahí dentro.
#
using Z_auxiliaresBRW, Z_auxiliaresGEN
#
using JLD
#
Γ = OS( "linux" ); # establecer SO a usar...es para mi lap si?
BRWname = "/home/isabel/Desktop/200819/control_01.brw"; # path file
# path para guardar los cachos
WORKpath = string( split( BRWname, "." )[ 1 ], Γ, "Process" ); 
checkpath( WORKpath ); # checar si el folder para guardar los cachos existe
#

#
ν = brw_things( BRWname ); # variables útiles
ξ = ν[ "Chs" ]; # canales en numero, coord1, coord2
Σ = ν[ "dset" ]; # path al dataset (NO READ! to much)
ζ = ν[ "Factor" ]; # variables para...
ο = ν[ "Offset" ]; # conversión de voltaje
Τ = Int( size( Σ, 1 )/size( ξ, 1 ) ); # verdadero tiempo registrado
#
ι = div_n_ab( Τ, 50, 100 ); # numero de cachos a cortar del exp. total
ω = Int( Τ/ι[ 1 ] ); # numero de frames finales (tamaño del cacho en frames)
ε = zeros( Int,ι[ 1 ], 2 ); # preallocate
ε[ :, 1 ] = collect( 1:ω:Τ ); # inicio y 
ε[ :, 2 ] = ε[ :, 1 ] .+ ω .- 1; # fin en frames de cada cacho (para cortar)
#

#
for i = 1:ι[ 1 ] # numero de Β a cortar (1-> 4096,1->ω)
    Β = zeros( size( ξ, 1 ), ω ); # preallocate
    # valores correspondientes al Β especifico
    # el canal 1,1 tiene el frame 1, 4097, 8193...etc
    β = collect( ( ε[ i, 1 ] - 1 ):1:ε[ i, 2 ] ); 
    for j = 1:ω
        # saca esos frames del machote seguido Σ, y ponlos en array en Β
        Β[ :, j ] = Σ[ ( β[ j ]*size( ξ, 1 ) ) + 1:( size( ξ, 1 )*β[ j + 1 ] ) ];
    end
    Β = ( ( Β.*ζ ) .+ ο ); # conversion a volaje
    ΒINname = string( WORKpath, Γ, "Β", add_zeros( 2, i ), ".jld" ); 
    save( ΒINname, "data", Β );
    println( add_zeros( 2, i ), " of ", ι[ 1 ], "...saved" );
end
#


