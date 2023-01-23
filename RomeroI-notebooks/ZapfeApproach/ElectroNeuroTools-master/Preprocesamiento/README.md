# Directorio de Preprocesamiento

Los archivos que se obtienen del aparato BioCAM de Brainwave tienen
la terminación *.brw. Estós son archivos hdf5, que contienen
una jerarquía de carpetas y archivos usando una notación tipo
diccionario para organizar los datos experimentales y otros
parámetros, tales como taza de muestreo (frequencia) y
escala de unidades (conversión a voltajes).

Los datos de diferentes formas de experimento requieren diferentes
preprocesamiento. Por ejemplo, los evocados requieren detectar
el "golpe" de estimulo y la actividad alrededor de este, y
promediar bajo varias instancias del estimulo, normalmente tres.
En cambio los experimentos espontaneos suelen ser más simples de
organizar.

Por otra parte, tal vez debido a actualizaciones de software del BioCAM,
no siempre están los datos en arreglos similares.
Hasta ahora hemos encontrado los casos en que los datos se encuentran
en un arreglo de 64 por 64 por <número de cuadros) y también en
uno de 4096 por <númmero de cuadros>.  Recomendamos una exploración interactiva
de los datos en un cuaderno Jupyter para tener una idea de como resulta el
caso.

