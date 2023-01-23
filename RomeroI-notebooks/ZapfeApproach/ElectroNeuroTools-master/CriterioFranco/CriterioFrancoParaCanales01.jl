#Implementemos el criterio (mas o menos complicado) de Franco
#Para encontrar canales con actividad neuronal


using MAT 

#Cualquier conjunto de datos sirve
DatosAparato=matopen("propagacion.mat")
matread("propagacion.mat")
LasMedicionesCrudas=read(DatosAparato, "RawMatrix")
#La parte que usamos para "jugar"
ParteInteresante=LasMedicionesCrudas[:,:,1300:2700];
medidasdatos=size(ParteInteresante)
tmax=medidasdatos[3]

datostest=reshape(ParteInteresante[9,30,:],tmax)

using FuncionFranco01

puchis=SigmaVentanas(datostest, 350, 35)

print(puchis)
