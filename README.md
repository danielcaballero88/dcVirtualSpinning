# dcVirtualSpinning
Pyhton package to produce virtual spinning microstructures

# Formatos de archivos de mallas

Los dos tipos de archivo (Mallacom y Mallasim) estan estructurados con etiquetas que empiezan con "\*".

## Mallacom 

Es una malla _as deposited_, el modulo mallacom permite obtener mallas de este tipo.
Luego, un programa externo de Fortran calcula las intersecciones y reescribe el archivo agregando "\_i" al nombre.

### Etiquetas

1. "\*Parametros"
    1. L: longitud dimensional de lado del RVE
    2. Dm: diametro medio de las fibras
    3. volfrac: fraccion de volumen de las fibras
    4. ls: longitud de segmento dimensional de las fibras
    5. devangmax: angulo de desviacion maximo entre segmentos (en grados)
    
2. "\*Coordenadas"
    1. n: numero de nodos en la malla
    2. (n lineas de):
        1. i: indice de nodo
        2. tipo: tipo de nodo (mas en seccion Nodos)
        3. x: coordenada x del nodo
        4. y: coordenada y del nodo
        
3. "\*Segmentos"
    1. n: numero de segmentos 
    2. (n lineas de):
        1. i: indice de segmento
        2. n0: primer nodo del segmento 
        3. n1: segundo nodo del segmento
        
4. "\*Fibras"
    1. n: numero de fibras
    2. (n lineas de):
        1. i: indice
        2. ls: longitud de segmentos
        3. d: diametro
        4. dth: angulo de desviacion maximo entre segmentos para esta fibra
        5. nsegs: numero de segmentos que componen la fibra
        6. con (integer(nsegs)): conectividad (indices de los segmentos)

5. "\*Capas"
    1. n
    2. (n lineas de):
        1. i
        2. nfibs
        3. con (integer(nfibs))
        
### Tipos de nodos

1. 0: Nodo interior. Son todos los nodos que surgen al depositar virtualmente segmentos lineales, entre los segmentos.
2. 1: Nodo frontera. Tambien ocurren de forma natural al depositar virtualmente las fibras, solamente que ocurren entre un segmento y la frontera.
3. 2: Nodo interseccion. Se calculan _a posteriori_ de depositar las fibras. Ocurren entre fibras de una misma capa o de capas adyacentes (por ahora, ya que esta modelizacion de las intersecciones puede cambiar).

## Mallasim

### Etiquetas

1. "\*Parametros"
    1. L
    2. Dm
    3. Nc
    4. nparam
    5. param (float(nparam))
    
2. 
