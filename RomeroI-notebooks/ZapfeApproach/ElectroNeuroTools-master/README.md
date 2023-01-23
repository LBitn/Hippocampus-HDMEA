Trabajos en Neurociencias desarrollados con el equipo del Dr. R. Gutiérrez
=====================================================================

-------------------------------------------------------------------

(c) Copyright  2015,2016,2017,2018,2019
Wilhelm Pablo Karel Zapfe Zaldivar.
Este software se distribuye con la licencia general pública GNU 3.0

--------------------------------------------------------------------


Aquí está un repositorio con los problemas tratados en el laboratorio 19
del CINVESTAV sede sur. El ataque a los problemas proviene de técnicas
numéricas basadas en razonamientos fíśicos clásicos y son extensiones
de técnicas tratadas anteriormente por otros científicos. Los separaremos en
carpetas dedicadas cada una a una parte independiente de los análisis
propuestos.

## CSDA

Técnicas de CSDA "clásico" o diferencial. Aquí podemos encontrar operadores
Laplacianos Numéricos para obtener el CSD de un LFP registrado con una
matriz de electrodos muy densa. La carpeta incluye técnicas de suavizado
de datos y de representación gráfica de los mismos. Los operadores
elegidos son los sugeridos por Lindenblat (operadores convexos).

## kCSDA

Está técnica para obtener el CSD es mucho más poderosa, basada en el trabajo de
Potworowski et al. Permite obtener el CSD de mallas irregulares de electrodos,
a través de funciones modelo ideales. El cálculo de los operadores
es la parte computacionalmente extensiva de esté trabajo.

## CenterOfMass

Aquí extendemos las ideas del cálculo de centro de masa de Elías y otros para
CSD. La aportación original consiste en separar la señal de pozos y fuentes
en componentes disjuntos, que permiten localizar la actividad en grupos
coherentes siempre cerca de los somas de las neuronas.


## CriterioFranco

Aquí probamos y extendemos varias ideas que Franco Ortiz propuso para
establecer correlaciones entre GD y CA.

## Graficando

Todas las rutinas para hacer figuras para paper, poster y presentaciones.

## Evocada

Rutinas para preprocesar (promediar) datos evocados, obtención de dCSD y CM.
Casi todas son extensiones sobre lo desarrollado en la carpeta CenterOfMass.
También hay rutinas que permiten separar electrodos saturados y electrodos con
actividad presuntamente fisiológica.

## PrimerosIntentos

Mayormente aquí encontramos un archivo muerto de código con los primeros
intentos de varios cálculos y cálculos tentativos de otras soluciones que aún
no hemos concluido con seriedad.


(c) Copyright  2015,2016,2017,2018,2019
Wilhelm Pablo Karel Zapfe Zaldivar.
Este software se distribuye con la licencia general pública GNU-3.0
