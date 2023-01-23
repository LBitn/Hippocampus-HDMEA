module MEATools

using HDF5, JLD, Dates, DelimitedFiles, Plots, HistogramThresholding, ImageContrastAdjustment,
      Suppressor, StatsBase, Statistics
export FindKeyRegex, FindKey, GetVarsHDF5, SearchDir, Zplot, Thresholding, ThresholdingPlots, FillingHolesCrux,
       FigureGroups, neighborgs, div_ab, SearchKey, Digital2Analogue, OneSegment,
DesatNegativePositive, ChunkSizeSpace, Donoho, ms2frames, SavePaths, GetGroups

# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
  SearchDir( path::String, key::String ) -> Vector{String}
      Find inside the given path, the files with the key word on their name and returns a vector of strings 
      with the full name of those files (complete path)
"""
SearchDir( path::String, key::String ) = filter( x -> endswith( x, key ), readdir( path; join = true ) );
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    FindKeyRegex( key::String, D::Dict{Any, Any} ) -> Vector{String}
        Searchs the entries of the Dictionary D for the key word using Word Match.
"""
function FindKeyRegex( key::String, D::Dict{Any, Any} )
    key = Regex( ".+$key.+" );
    ok = keys( D );
    aux = match.( key, ok );
    aux01 = aux[ aux .!= nothing ];
    okok = [ ];
    if !isempty( aux01 )
        for i in aux01
            push!( okok, i.match );
        end
        return okok
    else
        println( "There is no match with that key word into the Dicctionary" );
    end
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    FindKey( BRW::HDF5.File, key::String ) -> aux::String
        Find a key word into the HDF5 structure file and returns the path to the attribute or dataset
"""
function FindKey( BRW::HDF5.File, key::String )
    aux = [ ];
    for i in BRW
        if haskey( i, key )
            aux = string( i );
            groupaux = split( aux, " " )[ 2 ];
            aux = string( groupaux, "/", key );
        end
    end
    return aux
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    Z0( X::VecOrMat, nChs::Int64 ) -> Z::Matrix{Int64}
        using Plots
"""
function Z0( X::VecOrMat, nChs::Int64 )
    X = Int.( vec( X ) );
    Z = zeros( Int, nChs );
    n = Int( sqrt( nChs ) );
    Z[ X ] .= Z[ X ] .+ 1;
    Z = reverse( reshape( Z, n, n )', dims = 1 );
    return Int.( Z )
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    ZW( X::VecOrMat ) -> Z::Matrix{typeof(X)}
        using Plots
"""
function ZW( X::VecOrMat )
    X = vec( X );
    n = Int( sqrt( length( X ) ) );
    Z = reverse( reshape( X, n, n )', dims = 1 );
    return Z
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    Zplot( Z::Matrix, which::String, cm = :greys, nChs = 4096 ) -> F::Plot
        using Plots, MEATools.Z0, MEATools.ZW
"""
function Zplot( Z::VecOrMat, which::String, cm = :greys, nChs = 4096 )
    if which == "0"
        Z = Z0( Z, nChs ); 
    elseif which == "W"
        Z = ZW( Z );
    end
    F = heatmap( Z, aspect_ratio = 1, c = cm, axis = ( [ ], false ), wsize = ( 400, 400 ) );
    return F
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    GetGroups( BRW::HDF5.File, g::String ) -> GroupsN::Vector{String}
        Extract the groups form a BRW open file
"""
function GetGroupsHDF5( BRW::HDF5.File, g::String )
    GroupsN = [ ];
    try
        GroupsN = string.( g, "/", keys( BRW[ g ] ) );
    catch e
    end
    return GroupsN
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    GetAttr( BRW::HDF5.File, g::String ) -> AttrN::Vector{String}
        Extract the attributes form a BRW open file
        using HDF5
"""
function GetAttr( BRW::HDF5.File, g::String )
    AttrN = [ ];
    aux = attributes( BRW[ g ] );
    try
        AttrN = string.( g, "/", keys( aux ) );
    catch e
    end
    return AttrN
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    ExperimentSettings2Dict( Variables::Dict{Any, Any} ) -> Variables::Dict{Any, Any}
        Extracting the contents of the ExperimentSettings dictionary
"""
function ExperimentSettings2Dict( Variables::Dict{Any, Any} )
    ExperimentSettings = Variables[ "ExperimentSettings" ];
    t = split( ExperimentSettings, "\r\n" );
    t = replace.( t, "  " => "", "{" => "", "}" => "", '"' => "" );
    x = [ ];
    for i = 1:length( t )
        if !isempty( collect( eachmatch( r"[a-z]", t[ i ] ) ) )
            push!( x, true );
        else
            push!( x, false );
        end
    end
    t = t[ Bool.( x ) ]; t = split.( t, ": " );
    D = Dict( );
    for i in t
        if !( i[ 2 ] == "" )
            aux = i[ 2 ];
            try
                aux = replace( i[ 2 ], "," => " ", "[" => "", "[]" => "","]" => "", " " => "" )
            catch e
            end
            if ( aux != "" ) && ( aux != " " )
                aux = aux;
                try
                    aux = parse( Float64, aux );
                catch
                    aux = replace( aux, " " => "" );
                end
                D[ i[ 1 ] ] = aux;
            end
        end
    end
    delete!( Variables, "ExperimentSettings" );
    Variables = merge( Variables, D );
    return Variables
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    ExperimentDate2String( Variables::Dict{Any, Any} ) -> Variables::Dict{Any, Any}
        Extracting the date of the BRW creation
        using Dates
"""
function ExperimentDate2String( Variables::Dict{Any, Any} )
    X = Variables[ "ExperimentDateTime" ];
    Dt = split( X, ":" );
    Dt[ end ] = string( round( Int, parse( Float64, replace( Dt[ end ], r"[A-Z]" => "" ) ) ) );
    newDt = String( "" );
    for i in Dt
        newDt = string( newDt, ":", i );
    end
    newDt = newDt[ 2:end ];
    X = Dates.DateTime( newDt );
    Variables[ "ExperimentDateTime" ] = string( X );
    T = Dates.format( X, RFC1123Format );
    println( "Creation Date: ", T )
    return Variables
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    GetVarsHDF5( FILEBRW::String ) -> D::Dict{Any, Any}
        Extracts all contents from a HDF5 file of all versions, exept the Raw dataset, which only shows the 
        string with the location
        using HDF5, DelimitedFiles, JLD, MEATools.FindKeyRegex, MEATools.ExperimentSettings2Dict,
        MEATools.ExperimentDate2String, MEATools.GetAttr, MEATools.GetGroupsHDF5
"""
function GetVarsHDF5( FILEBRW::String )
    BRW = h5open( FILEBRW, "r" );
    Groups = keys( BRW );
    Groups = keys( BRW );
    AllGroups = [ ];
    AllAtributes = keys( attributes( BRW ) );
    while !isempty( Groups )
        GROUPS = [ ];
        ATTR = [ ];
        for g in Groups
            if typeof( BRW[ g ] ) == HDF5.Dataset
                push!( AllGroups, g )
            else
                push!( GROUPS, GetGroupsHDF5( BRW, g ) );
                AllAtributes = vcat( AllAtributes, GetAttr( BRW, g ) );
            end
        end
        Groups = vcat( GROUPS... );
        push!( AllGroups, Groups );
    end
    AllGroups = vcat( AllGroups... );
    Types = Vector{ String }( undef, length( AllGroups ) );
    for g = 1:length( AllGroups )
        Types[ g ] = string( typeof( BRW[ AllGroups[ g ] ] ) );
    end
    AllDataSets = AllGroups[ Types .== "HDF5.Dataset" ];
    aux = zeros( Int, length( AllDataSets ) );
    for i = 1:length( AllDataSets )
        aux[ i ] = length( BRW[ AllDataSets[ i ] ] );
    end
    Raw = AllDataSets[ aux .== maximum( aux ) ][ 1 ];
    NoRaw = AllDataSets[ aux .!= maximum( aux ) ];
    aux = aux[ aux .!= maximum( aux ) ];
    NoRaw = NoRaw[ aux .!= 0 ];
    D = Dict( );
    e = [ ];
    for g in NoRaw
        try 
            D[ g ] = Float64.( read( BRW[ g ] ) );
        catch e
            D[ g ] = read( BRW[ g ] );
        end
    end
    for g in keys( D )
        if length( D[ g ] ) .== 1
            D[ g ] = D[ g ][ 1 ];
        end
        try
            D[ g ] = Int64( D[ g ] );
        catch e
        end
    end
    for a in AllAtributes
        aux00 = basename( a );
        aux01 = dirname( a );
        if isempty( aux01 )
            D[ a ] = read_attribute( BRW, aux00 );
        else
            D[ a ] = read_attribute( BRW[ aux01 ], aux00 );
        end
    end
    D[ "Raw" ] = Raw;
    BRWsize = ( ( stat( FILEBRW ).size ) / 1000000 ) / 1024;
    D[ "BRWname" ] = BRW.filename;
    D[ "BRWsizeGB" ] = BRWsize;
    aux0 = FindKeyRegex( "Ch", D );
    aux1 = [ ];
    for i in aux0
        try
            push!( aux1, size( D[ i ], 1 ) );
        catch e
            push!( aux1, 0 );
        end
    end
    nChs = size( D[ aux0[ aux1 .== maximum( aux1 ) ][ ] ], 1 );
    D[ "nChs" ] = nChs;
    aux0 = FindKeyRegex( "Std", D );
    aux1 = [ ];
    for i in aux0
        try
            push!( aux1, size( D[ i ], 1 ) );
        catch e
            push!( aux1, 0 );
        end
    end
    STD = D[ aux0[ aux1 .== maximum( aux1 ) ][ ] ];
    D[ "STD" ] = STD;
    if ( "ExperimentSettings" in keys( D ) )
        D = ExperimentSettings2Dict( D );
        D = ExperimentDate2String( D );
    end
    X = String.( keys( D ) )[ values( D ) .== "null" ];
    for i in X
        delete!( D, i );
    end
    X = [ ];
    for i in keys( D )
        try
            if isempty( D[ i ] )
                push!( X, i );
            end
        catch e
        end
    end
    for i in X
        delete!( D, i );
    end
    x = [ ];
    for i in values( D )
        if length( i ) <= 100
            push!( x, 1 )
        else
            push!( x, 0 )
        end
    end
    NK = string.( keys( D ) )[ Bool.( x ) ];
    TXT = Dict( ); for i in NK; TXT[ i ] = D[ i ]; end
    PATHMAIN = split( FILEBRW, "." )[ 1 ]; mkpath( PATHMAIN );
    PATHINFO = joinpath( PATHMAIN, "Info" ); mkpath( PATHINFO );
    FILEVARS = joinpath( PATHINFO, "Variables.jld" );
    FILEVARSTXT = joinpath( PATHINFO, "Variables.txt" );
    FILEPATHS = joinpath( PATHINFO, "Paths.jld" );
    PATHBRWs = dirname( FILEBRW );
    PATHS = Dict(
        "PATHMAIN" => PATHMAIN,
        "PATHINFO" => PATHINFO,
        "PATHBRWs" => PATHBRWs
        );
    writedlm( FILEVARSTXT, TXT );
    save( FILEVARS, "Variables", D );
    save( FILEPATHS, "PATHS", PATHS );
    println( "You are now working on the new main path: ", PATHMAIN );
    println( "With the file: ")
    print( basename( BRW.filename ), " : ", D[ "Description" ] );
    println( "HDF5 file size: $BRWsize GB" );
    cd( PATHMAIN )
    close( BRW )
    return D, FILEPATHS
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
# using MeaTools.FindKeyRegex
function Digital2Analogue( Variables::Dict{ Any, Any }, DigitalValue::Matrix{UInt16} )
    SignalInversion = Variables[ FindKeyRegex( "SignalInversio", Variables )[ 1 ] ];
    MinVolt = Variables[ FindKeyRegex( "MinVol", Variables )[ 1 ] ];
    MaxVolt = Variables[ FindKeyRegex( "MaxVol", Variables )[ 1 ] ];
    BitDepth = Variables[ FindKeyRegex( "BitDept", Variables )[ 1 ] ];
    MVOffset = SignalInversion*MinVolt;
    ADCCountsToMV = ( SignalInversion * ( MaxVolt - MinVolt ) ) / ( 2^BitDepth );
    AnalogValue = @. MVOffset + ( DigitalValue * ADCCountsToMV )
    return AnalogValue
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function SearchKey( key, Variables )
    OKPATH = [ ];
    K = string.( keys( Variables ) );
    for k in K
        AUX = split( k, "/" );
        aux = match.( Regex(".+$key.+"), AUX );
        if !isempty( findall( aux .!= nothing ) )
             push!( OKPATH, join( AUX, "/" ) );
        end
    end
    if isempty( OKPATH )
        println( "There is no entry with the word $key" )
    else
        return OKPATH
    end
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    div_ab( n::Int64, lo::Int = 1, hi::Int = n )
        Divisors of the number n between the values "lo" and "hi", if they are not defined then it takes from 1
to n
"""
function div_ab( n::Int, lo::Int = 1, hi::Int = n )
    ρ = collect( 1:floor( Int, sqrt( n ) ) ) ; # the numbers one by one, from the square root
    σ1 = findall( n.%ρ .== 0 ); # square root divisors ( remainder = 0 )
    σ2 = Int.( ( n ) ./ ( σ1 ) ); # Take out the pairs (of 100, 2-50, 10-10, etc.)
    σ = sort( unique( vcat( σ1, σ2 ) ) ); # remove duplicates, concatenate, sort
    aux1 = @isdefined lo;
    aux2 = @isdefined hi;
    if aux1 && aux2
        rn = σ[ findall( hi .>= σ .>= lo ) ];
        if isempty( rn )
            println(" there is no divisors of $n between $lo and $hi" )
        else
            return rn
        end
    else
        return σ
    end
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    neighborgs( C::Int64, d::Int64 ) ->
        -> A = Array( ( d*2 ) + 1, ( d * 2 ) + 1 ), v = vec( 2*( ( d * 2 ) + 1 ) - 1 );
        The d-neighborhood is calculated from the channel (C) as a center
        A = array where C is the center and is in chip order
        v = same neighboring channels as A but in vector form and without C ( 8 channels )
"""
function neighborgs( C::Int64, d::Int64 )
    Layout = reverse( reshape( collect( 1:4096 ), 64, 64 )', dims = 1 );
    x_c = findall( Layout .== C )[ ][ 2 ]; y_c = findall( Layout .== C )[ ][ 1 ];
    aux = [ ( x_c - d ),( x_c + d ), ( y_c - d ), ( y_c + d ) ]
    aux[ aux .< 1 ] .= 1; aux[ aux .> 64 ] .= 64;
    A = Layout[ aux[ 3 ]:aux[ 4 ], aux[ 1 ]:aux[ 2 ] ];
    v = vec( A )[ vec( A ) .!= C ];
    return A, v
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function FigureGroups( grupos::Vector, loose::Vector = [ ], cm = :twilight )
    Z = zeros( Int, 4096 );
    for i = 1:size( grupos, 1 )
        Z[ grupos[ i ] ] .= floor( Int, log( length( grupos[ i ] ) ) ) + 2 ;
    end
    Z[ Int.( loose ) ] .= 1
    Z = reverse( reshape( Z, 64, 64 )', dims = 1 )
    F = heatmap(
        Z,
        aspect_ratio = 1,
        c = cm,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none );
    return F
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    FillingHolesCrux( Seleccion::Vector{Int64} ) -> NewSelection::Vector{Int64}
        using MeaTools.neighborgs
"""
function FillingHolesCrux( Seleccion::Vector{Int64} )
    Z = reverse( reshape( 1:4096, 64, 64 )', dims = 1 );
    arriba = Z[ 64, 2:( end - 1 ) ]; abajo = Z[ 1, 2:( end - 1 ) ];
    izquierda = Z[ 2:( end - 1 ), 1 ]; derecha = Z[ 2:( end - 1 ), 64 ];
    esquinas = vcat( Z[ 1, 1 ], Z[ 64, 64 ], Z[ 64, 1 ], Z[ 1, 64 ] );
    bordes = vcat( arriba, abajo, izquierda, derecha, esquinas );
    bordes_noesquinas = vcat( arriba, abajo, izquierda, derecha );
    P = true
    while P
        X1 = setdiff( setdiff( setdiff( 1:4096, Seleccion ) ), bordes );
        p1 = [ ];
        for i = 1:length( X1 )
            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
            CruzH = [ A[ 2, 1 ], A[ 2, 3 ] ];
            CruzV = [ A[ 1, 2 ], A[ 3, 2 ] ];
            x1 = length( CruzH[ CruzH .∈ [ Seleccion ] ] ) == 2;
            x2 = length( CruzV[ CruzV .∈ [ Seleccion ] ] ) == 2;
            if ( x1 || x2 )
                push!( p1, x )
            end
        end
        X1 = intersect( arriba, setdiff( setdiff( 1:4096, Seleccion ) ) );
        p2 = [ ];
        for i = 1:length( X1 )
            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
                Cruz = [ A[ 2, 1 ], A[ 2, 3 ] ];
                x1 = length( Cruz[ Cruz .∈ [ Seleccion ] ] ) == 2;
            if x1
                push!( p2, x )
            end
        end
        X1 = intersect( abajo, setdiff( setdiff( 1:4096, Seleccion ) ) );
        p3 = [ ];
        for i = 1:length( X1 )
            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
                Cruz = [ A[ 1, 1 ], A[ 1, 3 ] ];
                x1 = length( Cruz[ Cruz .∈ [ Seleccion ] ] ) == 2;
            if x1
                push!( p3, x )
            end
        end
        X1 = intersect( vcat( izquierda, derecha ), setdiff( setdiff( 1:4096, Seleccion ) ) );
        p4 = [ ];
        for i = 1:length( X1 )
            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
                Cruz = [ A[ 1, 2 ], A[ 3, 2 ] ];
                x1 = length( Cruz[ Cruz .∈ [ Seleccion ] ] ) == 2;
            if x1
                push!( p4, x )
            end
        end
        X1 = intersect( esquinas, setdiff( setdiff( 1:4096, Seleccion ) ) );
        p5 = [ ];
        for i = 1:length( X1 )
            x = X1[ i ];
            _, v = neighborgs( x, 1 );
            if length( v[ v .∈ [ Seleccion ] ] ) >= 2
               push!( p5, x )
            end
        end
        posibles = [ ];
        posibles = vcat( Int.( p1 ), Int.( p2 ), Int.( p3 ), Int.( p4 ), Int.( p5 ) );
        if isempty( posibles )
            P = false
        else
            Seleccion = vcat( posibles, Seleccion );
        end
    end
    return Int.( Seleccion )
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    ThresholdingPlots( Todos::Matrix{Int64} ) -> TF::Plot
    using MeaTools.Zplot, Plots
"""
function ThresholdingPlots( Todos::Matrix{Int64} )
    xmessage = string( "UnimodalRosin : ", length( findall( Todos[ :, 1 ] .== 1 ) )," channels" );
    T1 = Zplot( Todos[ :, 1 ], "W" ); T1 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "MinimumIntermodes: ", length( findall( Todos[ :, 2 ] .== 1 ) )," channels" );
    T2 = Zplot( Todos[ :, 2 ], "W" ); T2 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Intermodes: ", length( findall( Todos[ :, 3 ] .== 1 ) )," channels" );
    T3 = Zplot( Todos[ :, 3 ], "W" ); T3 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "MinimumError: ", length( findall( Todos[ :, 4 ] .== 1 ) )," channels" );
    T4 = Zplot( Todos[ :, 4 ], "W" ); T4 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Moments: ", length( findall( Todos[ :, 5 ] .== 1 ) )," channels" );
    T5 = Zplot( Todos[ :, 5 ], "W" ); T5 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Otsu: ", length( findall( Todos[ :, 6 ] .== 1 ) )," channels" );
    T6 = Zplot( Todos[ :, 6 ], "W" ); T6 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Entropy: ", length( findall( Todos[ :, 7 ] .== 1 ) )," channels" );
    T7 = Zplot( Todos[ :, 7 ], "W" ); T7 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Balanced: ", length( findall( Todos[ :, 8 ] .== 1 ) )," channels" );
    T8 = Zplot( Todos[ :, 8 ], "W" ); T8 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Yen: ", length( findall( Todos[ :, 9 ] .== 1 ) )," channels" );
    T9 = Zplot( Todos[ :, 9 ], "W" ); T9 = plot!( title = xmessage, cbar = :none );

    TF = plot(
        T1, T2, T3, T4, T5, T6, T7, T8, T9,
        layout = ( 3, 3 ), wsize = ( 700, 700 ), titlefont = ( 8, "arial" ) );
    return TF
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    Thresholding( W::VecOrMat ) -> Todos::Matrix{Int64}, T::Vector{Float64}
    using HistogramThresholding, ImageContrastAdjustment, Suppressor, StatsBase
"""
function Thresholding( W::VecOrMat )
    W = Float64.( vec( W ) );
    n = length( W );
    edges, conteo = HistogramThresholding.build_histogram( W, length( keys( countmap( W ) ) ) );
    t = zeros( 9 ); Todos = zeros( n, length( t ) );
    @suppress begin
        thr = find_threshold( UnimodalRosin( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 1 ] .= 1; t[ 1 ] = thr;
        thr = find_threshold( MinimumIntermodes( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 2 ] .= 1; t[ 2 ] = thr;
        thr = find_threshold( Intermodes( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 3 ] .= 1; t[ 3 ] = thr;
        thr = find_threshold( MinimumError( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 4 ] .= 1; t[ 4 ] = thr;
        thr = find_threshold( Moments( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 5 ] .= 1; t[ 5 ] = thr;
        thr = find_threshold( Otsu( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 6 ] .= 1; t[ 6 ] = thr;
        thr = find_threshold( Entropy( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 7 ] .= 1; t[ 7 ] = thr;
        thr = find_threshold( Balanced( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 8 ] .= 1; t[ 8 ] = thr;
        thr = find_threshold( Yen( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 9 ] .= 1; t[ 9 ] = thr;
    end
    return Int.( Todos ), t
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function ChunkSizeSpace( Variables::Dict{Any, Any}, limupper )
    NRecFrames = Variables[ SearchKey( "Frame", Variables )[ 1 ] ];
    σ = 5
    finalsize = Variables[ "BRWsizeGB" ] / σ;
    while ( finalsize > limupper )
        σ = σ + 1
        finalsize = Variables[ "BRWsizeGB" ] / σ;
    end
    while NRecFrames % σ != 0
        σ = σ + 1
    end
    println( "$σ segments of $finalsize GB each one" )
    return σ
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function DesatNegativePositive( limite::Float64, BIN::Matrix{Float64} )
    limite = size( BIN, 2 )*limite;
    AUX = findall( BIN .== minimum( BIN ) );
    AUX = getindex.( AUX, [ 1 2 ] );
    NegativeSaturatedChannels = unique( AUX[ :, 1 ] );
    NSFxC = Int.( values( countmap( AUX[ :, 1 ] ) ) );
    NSC = Int.( keys( countmap( AUX[ :, 1 ] ) ) );
    DiscardedChannelsN = NSC[ NSFxC .>= limite ];
    earth = NSC[ NSFxC .== maximum( NSFxC ) ][ 1 ];
    BIN[ earth, : ] .= 0;
    AUX = findall( BIN .== maximum( BIN ) );
    AUX = getindex.( AUX, [ 1 2 ] );
    PositiveSaturatedChannels = unique( AUX[ :, 1 ] );
    PSFxC = Int.( values( countmap( AUX[ :, 1 ] ) ) );
    PSC = Int.( keys( countmap( AUX[ :, 1 ] ) ) );
    DiscardedChannelsP = PSC[ PSFxC .>= limite ];
    DiscardedChannels = union( DiscardedChannelsN, DiscardedChannelsP );
    BIN[ DiscardedChannels, : ] .= 0;
    return earth, DiscardedChannels, BIN
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function OneSegment( Variables::Dict{Any, Any}, n::Int64, nSegments::Int )
    nChs = Variables[ "nChs" ];
    NRecFrames = Variables[ SearchKey( "Frame", Variables )[ 1 ] ];
    binlenght = Int( NRecFrames / nSegments );
    Raw = h5open( Variables[ "BRWname" ], "r" )[ Variables[ "Raw" ] ];
    BIN = Array{ UInt16 }( undef, nChs, binlenght );
    for frame = ( ( ( n - 1 ) * binlenght ) + 1 ): binlenght * n
        BIN[ :, ( frame - ( binlenght*( n - 1 ) ) ) ] .= Raw[ ( ( ( frame - 1 ) * nChs ) + 1 ): nChs * frame ];
    end
    return BIN
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
Donoho( x ) = ( median( abs.( x ) ) / 0.6745 );
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function ms2frames( time::Real, SamplingRate::Real )
    if time  != 0
        x = ceil( Int, ( time * SamplingRate ) / 1000 );
    else
        x = 1
    end
    return x
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function SavePaths( FILEPATHS::String )
    PATHS = load( FILEPATHS )[ "PATHS" ];
    vars = string.( names( Main )[ 4:end ] );
    pathsvars = String.( vars )[ findall( match.( r"PATH", String.( vars ) ) .!= nothing ) ];
    pathsfile = vcat( String.( keys( PATHS ) ), "FILEPATHS", "PATHS" );
    newpaths = pathsvars[ .!( pathsvars .∈ [ pathsfile ] ) ];
    if !isempty( newpaths )
        newdict = Dict( string( i ) => eval( Symbol( "$i" ) ) for i in newpaths );
        PATHS = merge( PATHS, newdict );
        save( FILEPATHS, "PATHS", PATHS );
        println( "Added ", keys( newdict ) );
    else
        println( "there is nothing new to add" )
    end
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function GetGroups( canales )
    A = neighborgs.( canales, 1 );
    A = getindex.( A, [ 1 ] );
    A = vec.( A );
    vecinos = [ ]
    for a in A
        push!( vecinos, a[ a .∈ [ canales ] ] );
    end
    loose = vcat( vecinos[ length.( vecinos ) .== 1 ]... );
    deleteat!( vecinos, findall( length.( vecinos ) .== 1 ) );
    NewGroups = Islas( vecinos );
    I = 1
    J = 2
    while I != J
        I = size( NewGroups, 1 )
        NewGroups = Islas( NewGroups );
        J = size( NewGroups, 1 );
    end
    return loose, NewGroups
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function Islas( vecinos )
    NewGroups = [ ];
    listaprueba = copy( vecinos );
    borrables = [ ];
    i = 1;
    while i <= size( vecinos, 1 )
        borrables = [ ];
        g = vecinos[ i ];
        for j in 2:size( vecinos, 1 )
            if !isempty( intersect( g, vecinos[ j ] ) )
                g = sort( union( g, vecinos[ j ] ) );
                push!( borrables, j )
            end
        end
        push!( borrables, i )
        push!( NewGroups, g )
        deleteat!( listaprueba, sort( unique( borrables ) ) )
        vecinos = copy( listaprueba );
    end
    return NewGroups
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
end # •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·• #


