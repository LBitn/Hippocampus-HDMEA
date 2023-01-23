#=
Usando el paquete ya listo de Julia, aprendamos a hacer FastICA primero
Valerí Makarov nos dice lo siguiente: Usa primero PCA, reduce la dimension, y haz ICA con eso. Lo más probable es que sea conveniente hacer todo en un archivo
que recuerde el modelo. Y luego reconstruir hacia atrás. De momento
Hagamos PCA aparte y luego ICA.
=#


using JLD, MultivariateStats

archivname=ARGS[1]
archivo=load(archivname)

# Variables generales


CompPrin=archivo["CompPrin"]

#=
Esto son solo los componentes principales, no la reconstruccion
Son mucho menos.
=#
(nchans,tmax)=size(CompPrin)



