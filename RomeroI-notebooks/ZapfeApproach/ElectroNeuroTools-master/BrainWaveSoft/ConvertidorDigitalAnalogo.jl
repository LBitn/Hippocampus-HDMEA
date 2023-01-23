module ConvertidorDigitalAnalogo


#= Esto convierte los valores en Uint de profundidad especifica a valores en milivoltios
usando como parametros los valores que el subarchivo /3BRecInfo/3BRecVars tiene guardados.

USO
Para llamar a la función, es necesario previamente cargar un archivo hdf5 usando por ejemplo:

using HDF5
DatosPrueba=h5open("/media/bodega/050815-hdf5/050815_1R0.brw", "r")

Despues de eso se puede llamar a la función, precedida del nombre del modulo:

test=ConvertidorDigitalAnalogo.ConvierteAMiliVolts(DatosPrueba)

entonces test tendra la misma forma que los datos contenidos en el archivo, pero en valores
de milivoltios medidos en el experimento.
 
=#


export ConvierteAMiliVolts

using HDF5

capacidad=500000


function ConvierteAMiliVolts(archivo::HDF5File)
    #Cargamos el pseudodirectorio donde estan los parametros del experimento
    ApuntaFactores=archivo["/3BRecInfo/3BRecVars"]
    Factores=read(ApuntaFactores)
    
    #vemos si el archivo es demasiado grande
    Frames=Factores["NRecFrames"][1]
    
    if Frames > capacidad
        println("Error: No puedo cargar todos los datos de un jalon")
        return 0;
    else
        #Cargamos Datos
        ApuntaDatos=archivo["/3BData/Raw"]
        #float es mas rapido que float16, float32 y float64 al parecer
        Aux=float(read(ApuntaDatos))
        #Aux=read(ApuntaDatos)       
        Signo=float(Factores["SignalInversion"][1])
        Divisor=float(2^Factores["BitDepth"][1])
        FactorConversion=(Factores["MaxVolt"][1]-Factores["MinVolt"][1])/(Divisor)*Signo
        Offset=float(Factores["MinVolt"][1])*Signo
        #println("Factor: ", FactorConversion, " Offset: ", Offset, "Signo: ", Signo)
        #Aux2=similar(Aux)
        Aux2=FactorConversion*Aux+Offset
        #Aux=0
        
    return Aux2
    end
end



end
