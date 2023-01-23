module Z_auxiliaresBRW
#-------------------#

using HDF5
function brw_things( file_brw::String ) # depende de brw file
    brw = h5open( file_brw, "r" );
    
    # channels and coords for each one
    Chs = read( brw[ "/3BRecInfo/3BMeaStreams/RawRanges" ] )[ "Chs" ]; 
    x = zeros( Int, size( Chs, 1 ), 3 );
    for i = 1:size( Chs, 1 )
        ## Los genios de 3Brain le volvieron a cambiar los nombres a las variables
        ## antes era data con dos entradas
        ## ahora es Row y Col
	##x[ i, : ] = [ i, Chs[ i ].data[ 1 ], Chs[ i ].data[ 2 ] ]; # Channels, [number coord1 coord2]
        x[ i, : ] = [ i, Chs[ i ].Row, Chs[ i ].Col ]; # Channels, [number coord1 coord2]
    end  
    
    ValidChs = read( brw[ "/3BData/3BInfo/3BNoise/ValidChs" ] ); # Valid Channels (software desition)
    y = zeros( Int, size( ValidChs, 1 ), 3 );
    for i = 1:size( ValidChs, 1 )
         ## Los genios de 3Brain le volvieron a cambiar los nombres a las variables
        ## antes era data con dos entradas
        ## ahora es Row y Col
	y[i,:] = [ ( ( ValidChs[ i ].Row - 1)*64 + ValidChs[ i ].Col ), ValidChs[ i ].Row, ValidChs[ i ].Col ];# Valid Channels, [number coord1 coord2]
    end
    
    RawEncodedTOC = read( brw[ "/3BData/RawEncodedTOC" ] )[ 2 ]; # numero real de frames registrados
    RecVars = read( brw[ "/3BRecInfo/3BRecVars" ] );
    Offset = RecVars[ "SignalInversion" ][ ]*RecVars[ "MinVolt" ][ ];
    Factor = RecVars[ "SignalInversion" ][ ]*( RecVars[ "MaxVolt" ][ ] - RecVars[ "MinVolt" ][ ] )/( 2^RecVars[ "BitDepth" ][ ] );
    NRecFrames = RecVars[ "NRecFrames" ][ ]; # numero propuesto de cuadros registrados
    SamplingRate = RecVars[ "SamplingRate" ][ ];
    dset = brw[ "3BData/RawEncoded" ];
    
    vars = Dict(
                "Offset"        => Offset,
                "Factor"        => Factor,
                #"NRecFrames"    => NRecFrames, # esta variable es inutil y erronea
                "SamplingRate"  => SamplingRate,
                "Chs"   	    => x,
		"ValidChs" 	    => y,
                "RawEncodedTOC" => RawEncodedTOC,
                "dset"          => dset,
                "MaxVolt"       => RecVars[ "MaxVolt" ][ ],
                "MinVolt"       => RecVars[ "MinVolt" ][ ]
        
    );

    return vars
    
end

#----#
export brw_things

#-#
end
