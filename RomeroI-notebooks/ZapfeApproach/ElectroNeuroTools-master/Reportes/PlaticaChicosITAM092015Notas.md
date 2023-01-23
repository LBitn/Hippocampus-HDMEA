De la medida eléctrica a lo subyacente.

El potencial eléctrico es un *efecto* debido a la distribución inhomogénea de
cargas eléctricas en el espacio. El potencial de una carga unitaria puntual
respecto a un cero colocado *infinitamente lejos* es $Q/r$.
Para distribuciones de cargas suaves, es posible considerar
la supèrposición continua de ellas como la siguiente integral, que
no es sino una suma formal sobre el espacio ocupado por las cargas.
El potencial de campo local resulta entonces una cantidad que decae como
1/r, mientras que la densidad de cargas mísmas decae como 1/r^3, lo
cual hace que sea mucho más densa cerca de las cargas que lo producen
y mucho menos importante lejos. Si pudieramos medir la densidad de cargas
instantanea o algo equivalente a ello, tendríamos algo mucho más preciso
en cuanto a la localización de la actividad neuronal.

#Del Potencial Local a la Densidad de Fuentes de Corriente.

Efectivamente existe una forma de invertir la expresión para obtener,
a partir del potencial, la densidad de cargas. El operador Laplaciano
en este caso efectivamente invierte la expresión. Un detalle que hay que considerar es que en este caso las cargas no son estáticas, sino que se desplazan,
y por razónes de interpretación fisiológica, más bien nos concentramos
en la densidad de fuentes de corrientes, que son proporcionales, por
conservación de la carga, al cambio de ésta en el tiempo. (ups, ahi va un
signo menos).

El problemilla con esto es que el Laplaciano es una idealización, como todo
operador diferencial. Nuestra escala mínima, es decir, nuestro epsilón, es
la escala de nuestras medidas en el mejor de los casos. Esto en realidad
tiene un lado bueno: las leyes diferenciales son aproximaciónes de grano grueso, y en realidad si pudieramos medir este límite nos topariamos con la estructura
celular en el camino, donde las cosas dejan de ser homogéneas e isotrópicas.
Así que en realidad estamos tomando límites pragmáticos, suficientemente chicos como para que la idealización matemática tenga sentido, y suficientemente grandes como para que pensemos que podemos promediar sobre volúmenes "infinitesimales" y olvidarnos de la estructura fina.

#Un Laplaciano Numérico

Muy bien, entonces tenemos que nuestro *"Laplaciano"* es una operación de
diferencias finitas alrededor del punto en que estamos interesados.
Dado que un verdadero Laplaciano es un operador invariante ante rotaciones,
debemos de poder obtener una representación numérica que sea lo más cercano
a esto, a pesar de que nuestra malla de datos es una estructura rectangular
que NO es invariante frente a rotaciones. Una buena sugerencia es
este operador convexo de Lindberg, cuya forma matricial es esta y su aplicación
es simplemente la convolución con la matriz de datos en cada instante
del tiempo. Esto es relativamente fácil de hacer operacionalmente y hace pocos
supuestos sobre la naturaleza de la distribución de cargas, excepto
isotropía e homogeneidad del medio externo (yo tengo mis dudas al respecto).

#Una aproximación inversa: iCSDA, kCSDA

Un método mucho más sofisticado es considerar lo que los físicos
llamamos el problema inverso, esto es, pensar que el voltaje que medimos
*ya es la solución* de un problema, pero nos falta plantear *bien* el problema.
Digamos que el Potencial de Campo es la solución de una distribución de cargas
medida también de forma discreta en N sitios. Entonces, lo que nos quedaría
por suponer es la distribución aproximada de las cargas en las dimensiones espaciales que nos faltan, y descubrir el kernel de una integral numérica.
Esta integral tiene que ver con el número de sitios que ya medimos.

#iCSDA, kCDA

Tenemos que resolver n^2 integrales que dependen cada una de n^2 sitios,
es decir, n^4 operaciones, y luego invertir esa matriz, que es una operación
de orden n^3 o en el mejor de los casos, n^2.3
Los investigadores que idearon este método para resolver el CSDA fueron muy lindos, y ya resolvieron el algorítmo para varias distribuciones posibles de
fuentes. Incluso pusieron a disposición su código e implementación en Matlab.
Desgraciadamente la mayor limitación de Matlab es su *overhead*. Si queremos
probar este método con la densidad de datos que nosotros tenemos, tenemos que reimplementer el método de forma astuta en un lenguaje más rápido y
usar todos los truquitos de A.L. para reducir las operaciones (simetría
de F, por ejemplo).

__cuello de botella__

# Un cero bien localizado.

¿Què otra virtud tiene la CSD sobre el LFP, aparte de una localización
precisa de los fenómenos? Pues tiene la virtud de que revela el signo
*real* de la carga. Entonces el *cero* que observamos en esta representación
separa dos grupos distintos con interpretación electríca clara:
uno es el atractor de cargas positivas, es decir, si dibujamos campos vectoriales, sería hacia donde se dirige la corriente, y se denomina el pozo. Los otros son las fuentes, lugares en donde se origina una sobrecarga positiva.
En esta representación el cero separa los dos conjuntos, que tienen por lo
tanto fronteras bien definidas. ¿qué podemos hacer con esto?

#Promedios vectoriales o centros de masa

La idea de poder seguir la trayectoria de actividad, o la sucesión de
centros de actividad se ha atacado de varias formas. Una de ellas fue
llevar a cabo promedios vectoriales de los centros de actividad, tomando
como ponderante las medidas del LFP sobre un arreglo de electrosos
al que consideramos un muestreo de un espacio vectorial. Esto es una aproximación válida, pero puede tener problemas de ambiguedad. Si queremos que la
suma de vectores sea convexa, la ponderación tiene que ser positiva.
En el LFP, lo que es positivo y negativo depende de constantes
arbitrarias. Lo que se llega ha considerar como masa, por ejemplo, es la
deflección hacia lo positivo en el LFP y se ignora la actividad que
causa deflección del LFP hacia el otro lado.

#Promedios vectoriales o centros de masa

En la representación CSD, tenemos dos ventajas: Uno, nuestro cero no es arbitrario, y podemos hacer la suma vectorial sobre cargas de uno u otro signo, siendo
la *densidad* de carga la ponderación adecuada, que garantiza una suma cóncava siempre. Si lo hicieramos directamente  podríamos obtener centros de actividad
un poco fuera de lugar, especialmente cuando los picos de la actividad están
separados. Las estructuras neuronales que estamos estudiando tienen sus
somas arreglados en figuras cóncavas, lo cual puede darnos "centros de actividad" sin ninguna validez fisiológica.

Peeero, el cero de nuestra función es una frontera clara. Así que si los picos
de la función son separados por elementos del otro conjunto,
cosa que típicamente ocurre), podemos separar el conjunto de todos los
pozos o fuentes en sus componentes disjuntos, y llevar a cabo el
promedio vectorial localizado *ahí* especificamente. De esta forma
el "centro" de actividad estaría contenido en el conjunto siempre
y sería una medida mucho más clara de un centro.
Esto por supuesto, requiere una rutina eficiente y una estructura de datos
flexible, cosas que yo no se implementar. Una vez resuelto ese detalle
podemos tener, punto por punto, el centro de masa de cada
*componente* de las densidades de carga y su magnitud.

__cuello de botella__

#Trayectorias de Actividad

¿Porqué insisto en la flexibilidad? Porque no nos interesa localizar
esos centros asi nomás. Para que podamos darles validez como representantes
de la actividad de las neuronas tenemos que poderlos seguir su movimioento.
Si estos aparentan tener cierta noción de continuidad, podríamos entonces
tal vez apreciar la actividad conjunta y secuencial de grupos de neuronas
juntos. Para ello debemos de poder distinguir entre un cm en un tiempo
dado y sus sucesores temporales efectivamente. Este es un problema computacional
complejo y la verdad, no tengo idea de como resolverlo efectivamente para
experimentos grandes.


__cuello de botella__




