# Global Imports
import csv
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
# Local Imports
from .Capas import Capas
from .Fibras import Fibras
from .Segmentos import Segmentos
from .Nodos import Nodos
from .Marco import Marco
from VirtualSpinning.aux import calcular_interseccion_entre_segmentos as calcular_interseccion
from VirtualSpinning.aux import find_string_in_file
from VirtualSpinning.aux import calcular_longitud_de_segmento
from VirtualSpinning.aux import calcular_angulo_de_segmento


MESH_PARAMS = ('L', 'D', 'vf', 'ls', 'dth', 'nc', 'fdo')
PI = np.pi


class Mallacom(object):
    def __init__(self, L, D, vf, ls, dth, nc, fdo=None):
        self.params = {
            'L': L,
            'D': D,
            'vf': vf,
            'ls': ls,
            'dth': dth,
            'fdo': fdo,
            'nc': nc
        }
        self.caps = Capas()  # lista vacia
        self.fibs = Fibras()  # lista vacia
        self.segs = Segmentos()  # lista vacia
        self.nods = Nodos()  # tiene dos listas vacias
        self.marco = Marco(L)

    def make_malla(self):
        """
        armo una malla
        """
        for _ic in range(self.params['nc']):
            self.make_capa()

    def make_capa(self, **kwargs):
        """
        armo una capa con fibras, todas van a armarse con los
        mismos parmetros dl y dtheta (se debe modificar para usar distribuciones)
        se depositan fibras hasta que se supera la fraccion de volumen dictada
        """

        # me fijo si reemplazo parametros globales de la malla
        cp = self.params.copy()  # cp = capa_params
        for key in kwargs: 
            assert key in MESH_PARAMS
            cp[key] = kwargs[key]

        # calculo el volumen de solido (ocupado por fibras) objetivo
        volc = cp['L']**2 * cp['D']  # volumen de la capa (fibras + vacio)
        vols_final = cp['vf'] * volc  # volumen de solido objetivo
        
        # --
        capa_con = []
        i = 0
        vols = 0.  # volumen de solido actual
        while vols < vols_final:
            i += 1
            j = self.make_fibra(**cp)
            if j == -1:
                i -= 1
            else:
                volf = self.calcular_volumen_de_una_fibra(j)
                vols += volf
                capa_con.append(j)

        # agrego la capa a la conectividad
        self.caps.add_capa(capa_con)
        
        # fin
        return 0

    def make_fibra(self, **params):
        """ tengo que armar una lista de segmentos
        nota: todos los indices (de nodos, segmentos y fibras)
        son globales en la malla, cada nodo nuevo tiene un indice +1 del anterior
        idem para segmentos y fibras
        los indices de los nodos, de los segmentos y de las fibras van por separado
        es decir que hay un nodo 1, un segmento 1 y una fibra 1
        pero no hay dos de misma especie que compartan indice """
        # ---
        # primero hago un segmento solo
        # para eso pongo un punto sobre la frontera del rve y el otro lo armo con un desplazamiento recto
        # tomo un angulo random entre 0 y pi, saliente del borde hacia adentro del rve
        # eso me da un nuevo segmento
        # agrego todas las conectividades
        # ---

        # Voy a ir guardando en una lista las coordenadas de los nodos
        coors = []

        # Variables locales mas comodas 
        fdo = params['fdo']
        ls = params['ls']
        dth = params['dth']
        L = params['L']
        D = params['D']

        # ----------
        # 1
        # Armo el primer segmento
        # primero busco un nodo en el contorno
        (x0, y0), b0 = self.marco.get_punto_random()
        if fdo is None:
            theta_fibra = np.random.rand() * PI
        else:
            theta_fibra = fdo()
        # Me fijo que la orientacion caiga en el rango periodico [0,pi)
        if theta_fibra == PI:
            theta_fibra = 0.
        elif theta_fibra > PI:
            raise ValueError("theta_fibra de una fibra no comprendido en [0,pi)")
        # Veo el cuadrante (-1: 0=hor, -2: pi/2=vert, 1: <pi/2, 2: >pi/2)
        cuad = self.get_cuadrante_theta_0_pi(theta_fibra)
        # Ahora me fijo la relacion entre cuadrante y borde
        check = self.check_if_fibra_alineada_a_borde(cuad, b0)
        if check == -1:
            return -1 # esta fibra no vale, esta alineada al borde
        else: 
            theta = self.get_theta_2pi(theta_fibra, cuad, b0)
        # Ya tengo el angulo del segmento
        dx = ls * np.cos(theta)
        dy = ls * np.sin(theta)
        coors.append([x0, y0])
        coors.append([x0 + dx, y0 + dy])
        # ----------

        # ----------
        # 2
        # Ahora agrego nuevos nodos (y segmenos) en un bucle
        # cada iteracion corresponde a depositar un nuevo segmento
        n = 1
        while True:
            # si el nodo anterior ha caido fuera del rve ya esta la fibra
            if self.check_fuera_del_RVE(coors[-1]):
                break
            n += 1
            # de lo contrario armo un nuevo segmento a partir del ultimo nodo
            # el angulo puede sufrir variacion
            theta = theta + dth * (2.0 * np.random.rand() - 1.0)
            # desplazamiento:
            dx = ls * np.cos(theta)
            dy = ls * np.sin(theta)
            # nuevo nodo
            x = coors[-1][0] + dx
            y = coors[-1][1] + dy
            coors.append([x, y])
        
        # ----------
        # 3
        # Aqui termine de obtener las coordenadas de los nodos que componen la fibra
        # si la fibra es muy corta la voy a descartar
        # para eso calculo su longitud de contorno
        # (luego la fibra se recorta asi que va a disminuir un poco)
        loco = ls * float(len(coors) - 1) 
        if loco < 0.3 * L:
            return -1
        # Voy a ensamblar la fibra como concatenacion de segmentos, que a su vez son concatenacion de dos nodos
        f_con = list()
        # agrego el primer nodo a la conectividad de nodos
        self.nods.add_nodo(coors[0], 1)
        for coor in coors[1:]:  # reocrro los nodos desde el nodo 1 (segundo nodo)
            self.nods.add_nodo(coor, 0)
            nnods = len(self.nods)
            s0 = [nnods - 2, nnods - 1]
            self.segs.add_segmento(s0, self.nods.r)
            nsegs = len(self.segs)
            f_con.append(nsegs - 1)
        # al final recorto la fibra y la almaceno
        self.nods.tipos[-1] = 1
        self.fibs.add_fibra(f_con, ls, D, dth, loco)
        # trim la fibra recientemente agregada (en la numeracion es la ultima)
        self.trim_fibra_at_frontera(len(self.fibs.con) - 1)  
        # ----------

        # fin
        return len(self.fibs.con) - 1  # devuelvo el indice de la fibra

    def calcular_loco_de_una_fibra(self, f):
        """ calcula la longitud de contorno de una fibra """
        loco = 0.
        for seg in self.fibs.con[f]:
            n0, n1 = self.segs.con[seg]
            r0 = self.nods.r[n0]
            r1 = self.nods.r[n1]
            lseg = calcular_longitud_de_segmento(r0, r1)
            loco += lseg
        return loco

    def calcular_volumen_de_una_fibra(self, f):
        """ calcula el volumen ocupado por una fibra """
        loco = self.calcular_loco_de_una_fibra(f)
        D = self.fibs.D[f]
        return loco * PI * D * D / 4.

    def calcular_fraccion_de_volumen_de_una_capa(self, capcon):
        """ calcula la fraccion de volumen de una capa como
        el volumen ocupado por fibras
        dividido el volumen total de la capa """
        volfs = 0
        for f in capcon:  # recorro las fibras de la capa
            volf = self.calcular_volumen_de_una_fibra(f)
            volfs += volf
        # el volumen total de la capa es:
        # TODO: el calculo deberia hacerse con los parametros de la capa
        # pero esos no estan disponibles ahora mismo, deberia guardarlos en la capa
        volc = self.params['L']**2 * self.params['D']
        # luego la fraccion de volumen
        fracvol = volfs / volc
        return fracvol

    def calcular_orientacion_extremo_extremo_de_una_fibra(self, f):
        fcon = self.fibs.con[f]
        s0 = fcon[0]
        s1 = fcon[1]
        n0 = self.segs.con[s0][0]
        n1 = self.segs.con[s1][1]
        r0 = self.nods.r[n0]
        r1 = self.nods.r[n1]
        theta_2pi = calcular_angulo_de_segmento(r0, r1)
        if theta_2pi >= PI:
            theta = theta_2pi - PI
        else:
            theta = theta_2pi
        if theta_2pi < 0.:
            raise ValueError
        return theta

    @staticmethod
    def get_cuadrante_theta_0_pi(theta):
        if theta < PI * 1.0e-8:
            return -1  # direccion horizontal
        elif np.abs(theta - PI * 0.5) < 1.0e-8:
            return -2  # direccion vertical
        elif theta < PI * 0.5:
            return 1  # primer cuadrante
        else:
            return 2  # segundo cuadrante

    @staticmethod
    def check_if_fibra_alineada_a_borde(cuad, b0):
        if cuad == -1:  # fibra horizontal
            if b0 in (0, 2):
                return -1  # esta fibra no vale, es horizontal sobre un borde horizontal
        elif cuad == -2:  # fibra vertical
            if b0 in (1, 3):
                return -1 # fibra vertical sobre borde vertical, no va
        return 0

    @staticmethod
    def get_theta_2pi(th_0pi, cuad, b0):
        """ a partir del angulo de fibra comprendido en [0, pi)
        y el borde del cual parte la fibra, se calcula el angulo
        que corresponde entre [0, 2pi) 
        Parametros:
            th_0pi: angulo de fibra en [0, pi) 
            cuad: cuadrante de la fibra: -1, -2, 1, 2 = horizontal, vertical, primer, segundo
            b0: integer con el borde, 0,1,2,3 = bottom, right, top, left
        Returns: 
            theta: angulo del versor de la fibra en [0, 2pi) hacia dentro del rve
        """ 
        if cuad == -1:  # fibra horizontal
            if b0 in (0, 2):
                raise Exception('Fibra horizontal alineada con borde!') 
            elif b0 == 1:
                theta = PI
            else:  # b0 == 3
                theta = 0.
        elif cuad == -2:  # fibra vertical
            if b0 in (1, 3):
                raise Exception('Fibra vertical alineada con borde!')
            elif b0 == 0:
                theta = 0.5 * PI
            else:  # b0 == 2
                theta = 1.5 * PI
        elif cuad == 1:  # primer cuadrante
            if b0 in (0, 3):
                theta = th_0pi
            else:  # b0 in(1,2)
                theta = th_0pi + PI
        else:  # cuad == 2 segundo cuadrante
            if b0 in (0, 1):
                theta = th_0pi
            else:  # b0 in (2,3)
                theta = th_0pi + PI
        return theta

    def check_fuera_del_RVE(self, r):
        x = r[0]
        y = r[1]
        L = self.params['L']
        if x <= 0 or x >= L or y <= 0 or y >= L:
            return True
        else:
            return False

    def trim_fibra_at_frontera(self, j):
        """ subrutina para cortar la fibra que ha salido del rve """
        # debo cortar la ultima fibra en su interseccion por el rve
        # para eso calculo las intersecciones de los nodos con los bordes
        # coordenadas del ultimo segmento de la fibra de conectividad fib_con
        fib_con = self.fibs.con[j]
        s = fib_con[-1]
        rs0, rs1 = self.nods.r[self.segs.con[s]]  # coordenadas xy de los nodos del segmento s
        # pruebo con cada borde
        for b in range(4):  # recorro los 4 bordes
            # puntos del borde en cuestion
            rb0, rb1 = self.marco.get_side_nodes(b) # coordenadas xy de los dos nodos del borde b
            interseccion = calcular_interseccion(rs0, rs1, rb0, rb1)
            if interseccion is None:  # no hubo interseccion
                continue  # con este borde no hay interseccion, paso al que sigue
            else:  # hubo interseccion
                in_r, in_tipo, *_ = interseccion
                if in_tipo == 2:  # interseccion en el medio
                    try:  # tengo que mover el ultimo nodo y por lo tanto cambia el segmento
                        self.segs.mover_nodo(s, 1, self.nods.r, in_r)
                        rs1 = in_r
                    except ValueError as err:
                        # TODO: limpiar esto por el amor de la virgen y todos los santos
                        print("Error en trim_fibra_at_frontera")
                        print(err)
                        print(fib_con, b, interseccion)
                        quit()
                else:  # interseccion coincide con uno o dos extremos
                    # en este caso solo me importa el segundo nodo del primer segmento (seg 0)
                    # porque el segmento 1 es el borde, y el primer nodo del seg 0 siempre deberia estar dentro del rve
                    # (o en el borde a lo sumo si se trata de una fibra de un solo segmento)
                    # y en ese caso no hay nada que hacer! puesto que el nodo ya esta en el borde
                    # TODO: programar este caso, no es complicado
                    pass

    def mover_nodo(self, n, s, f, new_r):
        """
        muevo un nodo (cambio su posicion), para eso tambien doy 
        a que segmentos y a que fibras (puede ser mas de 1?) pertenece

        Args:
            n: integer, indice del nodo
            s: lista de int con indices de segmentos
            f: lista de int con indices de fibras

        Returns:
            nada, solo muta a los nodos, segmentos y fibras
        """
        self.nods.r[n] = new_r 
        old_longs = []
        new_longs = []
        for js in s:
            # para cada segmento varia su long y theta (este metodo lo arregla)
            self.segs.actualizar_segmento(js, self.nods.r)
        for jf in f: 
            # para la fibra solo cambia su longitud de contorno (loco)


    def cambiar_capas(self, new_ncapas):
        """ un mapeo de las fibras en un numero de capas diferente """
        # me fijo cuantas fibras van a entrar en cada capa
        nfibras = len(self.fibs.con)
        nf_x_capa = int(nfibras / new_ncapas)
        # armo una nueva conectividad de capas
        capas_con = list()
        for c in range(new_ncapas - 1):  # -1 porque la ultima capa la hare aparte
            print(c)
            capa_con = range(c * nf_x_capa, (c + 1) * nf_x_capa)
            capas_con.append(capa_con)
        # la ultima capa puede tener alguna fibra de mas
        c = new_ncapas - 1
        capa_con = range(c * nf_x_capa, (c + 1) * nf_x_capa)
        capas_con.append(capa_con)
        self.caps.set_capas_listoflists(capas_con)

    def calcular_conectividad_de_interfibras(self):
        """ 
        Funcion para calcular las interfibras (fibras simplificadas)
        Se recorren las fibras (sus conectividades) y se arman las 
        interfibras (sus conectividades) descartando los nodos tipo 0
        y cambiando de interfibra en los nodos tipo 2
        -
        ojo son diferentes a las subfibras de una malla simplificada
        aqui las interfibras son concatenaciones de segmentos entre nodos interseccion
        en una ms las subfibras son una simplificacion de una fibra 
        dando solamente los nodos extremos y enrulamiento """
        infbs_con = list()  # conectividad: lista de listas de segmentos
        for fcon in self.fibs.con:  # recorro las fibras
            # cada fibra que empieza implica una nueva interfibra
            infb = list()  # conectividad de la interfibra: lista de segmentos
            # tengo que ir agregando segmentos hasta toparme con un nodo interseccion o frontera
            for s in fcon:  # recorro los segmentos de la fibra f
                scon = self.segs.con[s]
                _n0, n1 = scon
                # agrego el segmento s a la interfibra
                infb.append(s)
                # si el ultimo nodo de s es interseccion o frontera aqui termina la interfibra
                if self.nods.tipos[n1] in (1, 2):
                    infbs_con.append(infb)  # agrego la interfibra a la conectividad
                    infb = list()  # preparo una nueva interfibra vacia para continuar agregando segmentos
        # aqui ya deberia tener una conectividad terminada
        return infbs_con

    def calcular_orientaciones(self):
        """ calcular las orientaciones de las fibras de toda la malla """
        thetas_f = list()
        l = len(self.fibs.con)
        for f in range(l):
            theta_f = self.calcular_orientacion_extremo_extremo_de_una_fibra(f)
            thetas_f.append(theta_f)
        return thetas_f

    def calcular_distribucion_de_orientaciones(self, rec_orientaciones=None, n=10):
        """ calcula la distribucion de orientaciones en la malla
        contando las frecuencias en los bins """
        # obtengo las orientaciones de todas las fibras
        if rec_orientaciones is None:
            phis = self.calcular_orientaciones()
        else:
            phis = rec_orientaciones
        phis = np.array(phis, dtype=float)
        #
        conteo, x_edges = np.histogram(phis, bins=n, range=(0., PI))
        delta = (PI - 0.) / float(n)
        # pdf = conteo / float(np.sum(conteo)) / delta
        x = x_edges[:-1] + 0.5 * delta
        return x, delta, conteo

    def calcular_enrulamientos(self):
        """ calcular para todas las fibras sus longitudes de contorno y
        sus longitudes extremo a extremos (loco y lete)
        y calcula el enrulamiento como lamr=loco/lete """
        lamsr = []
        for fcon in self.fibs.con:  # recorro las fibras del rve
            loco = 0.
            for s in fcon:  # recorro los segmentos de cada fibra
                scon = self.segs.con[s]
                n0, n1 = scon
                r0 = self.nods.r[n0]
                r1 = self.nods.r[n1]
                try:
                    loco += calcular_longitud_de_segmento(r0, r1)
                except ValueError:
                    raise ValueError("Error, segmento de longitud nula!!")
            n_ini = self.segs.con[fcon[0]][0]
            n_fin = self.segs.con[fcon[-1]][1]
            r_ini = self.nods.r[n_ini]
            r_fin = self.nods.r[n_fin]
            try:
                lete = calcular_longitud_de_segmento(r_ini, r_fin)
            except ValueError:
                raise ValueError("Error, lete de longitud nula!!")
            lamsr.append(loco / lete)
        return lamsr

    def calcular_distribucion_de_enrulamiento(self, rec_lamsr=None, lamr_min=None, lamr_max=None, n=10, binwidth=None):
        """ calcular la distribucion de enrulamientos (pdf)
        para eso calculo el histograma y luego normalizo con el area
        parametros de entrada
        rec_lamsr: valores de lamsr para todas las fibras (por si ya lo tengo asi no lo calculo al dope de nuevo)
        lamr_min: valor minimo de lamr para hacer el histograma (ojo pueden quedar fibras afuera)
        lamr_max: valor maximo de lamr para hacer el histograma (ojo pueden quedar fibras afuera)
        n: numero de bins
        parametros de salida: x, delta, pdf
        x: valores medios de cada bin
        delta: ancho de cada bin (son todos iguales)
        pdf: valor de pdf
        """
        if rec_lamsr is None:
            lamsr = self.calcular_enrulamientos()
        else:
            lamsr = rec_lamsr
        lamsr = np.array(lamsr, dtype=float)
        if lamr_min is None: lamr_min = np.min(lamsr)
        if lamr_max is None: lamr_max = np.max(lamsr)
        # me fijo si hay imposicion de ancho de bin, 
        # entonces calculo con eso el numero de bins (es el entero mas cercano)
        if binwidth is not None: n = int((lamr_max - lamr_min) / binwidth + 0.5) 
        # calculo el histograma
        conteo, x_edges = np.histogram(lamsr, bins=n, range=(lamr_min, lamr_max))
        delta = (lamr_max - lamr_min) / n
        # pdf = conteo / float(np.sum(conteo)) / delta
        x = x_edges[:-1] + 0.5 * delta
        return x, delta, conteo

    def get_histograma_lamr(self, 
                            lamr_min=None, lamr_max=None, 
                            nbins=5, binwidth=None, opcion="fibras",
                            csv_file=False):
        if opcion == "fibras":
            lamsr = self.calcular_enrulamientos()
        elif opcion == "interfibras":
            lamsr = self.calcular_enrulamientos_de_interfibras()
        else:
            raise ValueError
        #
        lrs, dlr, conteo = self.calcular_distribucion_de_enrulamiento(rec_lamsr=lamsr, lamr_min=lamr_min,
                                                                      lamr_max=lamr_max, n=nbins, binwidth=binwidth)
        frecs = np.array(conteo, dtype=float) / float(np.sum(conteo))
        pdf = frecs / dlr
        # si hay opcion de guardar en archivo csv el histograma, lo hago:
        if csv_file:
            header = ['lamrs', 'dlamr', 'count', 'frec', 'pdf']
            columns = [lrs, [dlr] * len(lrs), conteo, frecs, pdf]
            rows = [list(row_tuple) for row_tuple in zip(*columns)]
            with open(csv_file, 'w') as file:
                writer = csv.writer(file)
                writer.writerow(header)
                writer.writerows(rows)
        return lrs, dlr, conteo, frecs, pdf

    def get_histograma_orientaciones(self, nbins=5, opcion="fibras", csv_file=False):
        if opcion == "fibras":
            rec_thetas = self.calcular_orientaciones()  # todos los angulos de las fibras
            # thetas a continuacion son los angulos medio de cada bin
        elif opcion == "interfibras":
            raise NotImplementedError
        else:
            raise ValueError
        #
        thetas, dth, conteo = self.calcular_distribucion_de_orientaciones(rec_orientaciones=rec_thetas, n=nbins)
        frecs = np.array(conteo, dtype=float) / float(np.sum(conteo))
        pdf = frecs / dth
        # si hay opcion de guardar en archivo csv el histograma, lo hago:
        if csv_file:
            header = ['theta', 'dth', 'count', 'frec', 'pdf']
            columns = [thetas, [dth] * len(thetas), conteo, frecs, pdf]
            rows = [list(row_tuple) for row_tuple in zip(*columns)]
            with open(csv_file, 'w') as file:
                writer = csv.writer(file)
                writer.writerow(header)
                writer.writerows(rows)
        return thetas, dth, conteo, frecs, pdf

    def calcular_enrulamientos_de_interfibras(self, rec_infbs_con=None):
        """
        Calcular para todas las interfibras sus longitudes de contorno y
        sus longitudes extremo a extremos (loco y lete)
        y calcula el enrulamiento como lamr=loco/lete
        
        Args:
            rec_infbs_con: recorded interfibers conectivity, se utiliza para 

        Returns:
            lamr: lista con los valores de enrulamiento para las interfibras
        """
        
        # inicio una lista vacia
        lamsr = []

        # si no ingrese la conectividad, la calculo (ineficiente)
        if rec_infbs_con is None:
            infbs_con = self.calcular_conectividad_de_interfibras()
        else: 
            infbs_con = rec_infbs_con

        # recorro las interfibras (fibras interectadas) del rve
        for infb_con in infbs_con: 
            loco = 0.
            for s in infb_con:  # recorro los segmentos de cada interfibra
                scon = self.segs.con[s]
                n0, n1 = scon
                r0 = self.nods.r[n0]
                r1 = self.nods.r[n1]
                loco += calcular_longitud_de_segmento(r0, r1)
            n_ini = self.segs.con[infb_con[0]][0]
            n_fin = self.segs.con[infb_con[-1]][1]
            r_ini = self.nods.r[n_ini]
            r_fin = self.nods.r[n_fin]
            lete = calcular_longitud_de_segmento(r_ini, r_fin)
            lamsr.append(loco / lete)

        # return
        return lamsr

    def calcular_distribucion_de_enrulamiento_de_interfibras(self, rec_lamsr=None, lamr_min=None, lamr_max=None, n=10):
        """ calcular la distribucion de enrulamientos (pdf)
        para eso calculo el histograma y luego normalizo con el area
        parametros de entrada
        rec_lamsr: valores de lamsr para todas las fibras (por si ya lo tengo asi no lo calculo al dope de nuevo)
        lamr_min: valor minimo de lamr para hacer el histograma (ojo pueden quedar fibras afuera)
        lamr_max: valor maximo de lamr para hacer el histograma (ojo pueden quedar fibras afuera)
        n: numero de bins
        parametros de salida: x, delta, pdf
        x: valores medios de cada bin
        delta: ancho de cada bin (son todos iguales)
        pdf: valor de pdf
        """
        if rec_lamsr is None:
            lamsr = self.calcular_enrulamientos_de_interfibras()
        else:
            lamsr = rec_lamsr
        lamsr = np.array(lamsr, dtype=float)
        if lamr_min is None:
            lamr_min = np.min(lamsr)
        if lamr_max is None:
            lamr_max = np.max(lamsr)
        delta = (lamr_max - lamr_min) / n
        conteo, x_edges = np.histogram(lamsr, bins=n, range=(lamr_min, lamr_max))
        delta = (lamr_max - lamr_min) / n
        pdf = conteo / float(np.sum(conteo)) / delta
        x = x_edges[:-1] + 0.5 * delta
        return x, delta, pdf

    def guardar_en_archivo(self, archivo="Malla.txt"):
        fid = open(archivo, "w")
        # ---
        # primero escribo L, ls y dtheta
        fid.write("*Parametros (L, Dm, volfrac, ls, devangmax) \n")
        fid.write("{:20.8f}\n".format(self.params['L']))
        fid.write("{:20.8f}\n".format(self.params['D']))
        fid.write("{:20.8f}\n".format(self.params['vf']))
        fid.write("{:20.8f}\n".format(self.params['ls']))
        fid.write("{:20.8f}\n".format(self.params['dth'] * 180. / PI))  # lo escribo en grados
        # ---
        # escribo los nodos: indice, tipo, y coordenadas
        dString = "*Coordenadas \n" + str(len(self.nods.r)) + "\n"
        fid.write(dString)
        for n in range(len(self.nods.r)):
            dString = "{:12d}".format(n)
            dString += "{:2d}".format(self.nods.tipos[n])
            dString += "".join("{:+17.8e}".format(val) for val in self.nods.r[n]) + "\n"
            fid.write(dString)
        # ---
        # sigo con los segmentos: indice, nodo inicial y nodo final
        dString = "*Segmentos \n" + str(len(self.segs.con)) + "\n"
        fid.write(dString)
        for s in range(len(self.segs.con)):
            n0, n1 = self.segs.con[s]  # conectividad del segmento (nodo inicial n0 y nodo final n1)
            fmt = "{:12d}" * 3
            dString = fmt.format(s, n0, n1) + "\n"
            fid.write(dString)
        # ---
        # sigo con las fibras: indice, ls, d, dtheta, nsegs_f, y segmentos (conectividad)
        dString = "*Fibras \n" + str(len(self.fibs.con)) + "\n"
        fid.write(dString)
        for f, fcon in enumerate(self.fibs.con):
            dString = "{:12d}".format(f)  # indice
            dString += "{:17.8e}{:17.8e}{:17.8e}".format(self.fibs.dl[f], self.fibs.D[f],
                                                         self.fibs.dth[f])  # ls, d y dtheta
            dString += "{:12d}".format(len(fcon))  # numero de segmentos en la fibra
            dString += "".join("{:12d}".format(val) for val in fcon) + "\n"  # conectividad
            fid.write(dString)
        # termino con las capas: indice y fibras (conectividad):
        dString = "*Capas \n" + str(len(self.caps.con)) + "\n"
        fid.write(dString)
        for c, ccon in enumerate(self.caps.con):
            dString = "{:12d}".format(c)  # indice
            dString += "{:12d}".format(len(ccon))  # numero de fibras en la capa
            dString += "".join("{:12d}".format(val) for val in ccon) + "\n"  # conectividad
            fid.write(dString)
        # ---
        # termine
        fid.close()

    @classmethod
    def leer_de_archivo(cls, archivo="Malla.txt"):
        fid = open(archivo, "r")
        # primero leo los parametros
        target = "*parametros"
        _ierr = find_string_in_file(fid, target, True)
        L = float(next(fid))
        Dm = float(next(fid))
        volfrac = float(next(fid))
        if volfrac > 1.:
            volfrac = int(volfrac + .5)
        ls = float(next(fid))
        devangmax = float(next(fid))  # en grados
        devangmax = devangmax * PI / 180.
        # luego busco coordenadas
        target = "*coordenadas"
        _ierr = find_string_in_file(fid, target, True)
        num_r = int(next(fid))
        coors = []
        tipos = []
        for i in range(num_r):
            _j, t, x, y = (float(val) for val in next(fid).split())
            tipos.append(int(t))
            coors.append([x, y])
        # luego los segmentos
        target = "*segmentos"
        _ierr = find_string_in_file(fid, target, True)
        num_s = int(next(fid))
        segs = []
        for i in range(num_s):
            _j, n0, n1 = (int(val) for val in next(fid).split())
            segs.append([n0, n1])
        # luego las fibras
        target = "*fibras"
        _ierr = find_string_in_file(fid, target, True)
        num_f = int(next(fid))
        fibs = []
        lss = []
        ds = []
        dthetas = []
        for i in range(num_f):
            svals = next(fid).split()
            _j = int(svals[0])
            ls = float(svals[1])
            d = float(svals[2])
            dtheta = float(svals[3])
            # _nsegsf = int(svals[4])  # unused
            fcon = [int(val) for val in svals[5:]]
            fibs.append(fcon)
            lss.append(ls)
            ds.append(d)
            dthetas.append(dtheta)
        # luego la capas
        target = "*capas"
        _ierr = find_string_in_file(fid, target, True)
        num_c = int(next(fid))
        caps = []
        for c in range(num_c):
            svals = next(fid).split()
            _j = int(svals[0])
            # _nfibsc = int(svals[1])  # unused
            ccon = [int(val) for val in svals[2:]]
            caps.append(ccon)
        # ahora que tengo todo armo el objeto
        params = {
            'L': L,
            'D': Dm, 
            'vf': volfrac, 
            'ls': ls,
            'dth': devangmax,
            'nc': num_c
        }
        malla = cls(**params)
        # le asigno los nodos
        for i in range(num_r):
            malla.nods.add_nodo(coors[i], tipos[i])
        # le asigno los segmentos
        for i in range(num_s):
            s_con = segs[i]
            try:
                malla.segs.add_segmento(s_con, coors)
            except ValueError:
                print("Segmento de long nula")
                # raise ValueError("Error, segmento de longitud nula!!")
        # le asigno las fibras
        for i in range(num_f):
            f_con = fibs[i]
            ls = lss[i]
            d = ds[i]
            dtheta = dthetas[i]
            malla.fibs.add_fibra(f_con, ls, d, dtheta)
        # le asigno las capas
        for c in range(num_c):
            c_con = caps[c]
            malla.caps.add_capa(c_con)
        # listo
        return malla

    def pre_graficar_capas(self, fig, ax, byn=True):
        nc = len(self.caps.con)
        if byn:
            mi_colormap = cm.get_cmap('gray')
            # mi_colormap = plt.cm.gray
        else:
            # mi_colormap = plt.cm.rainbow
            mi_colormap = cm.get_cmap('rainbow')
        sm = plt.cm.ScalarMappable(cmap=mi_colormap, norm=plt.Normalize(vmin=0, vmax=nc - 1))
        # dibujo las fibras (los segmentos)
        # preparo las listas, una lista para cada fibra
        xx = [[] for f in self.fibs.con]
        yy = [[] for f in self.fibs.con]
        grafs = []
        for c, c_con in enumerate(self.caps.con):  # recorro las capas
            for f in c_con:  # recorro las fibra de la capa
                f_con = self.fibs.con[f]
                # antes de recorrer los segmentos de cada fibra
                # el primer nodo del primer segmento lo agrego antes del bucle
                s = f_con[0]  # primer segmento de la fibra f
                n = self.segs.con[s][0]  # primer nodo del segmento s
                r = self.nods.r[n]  # coordenadas de ese nodo
                xx[f].append(r[0])
                yy[f].append(r[1])
                for s in f_con:  # recorro los segmentos de la fibra f
                    s_con = self.segs.con[s]
                    n = s_con[1]  # ultimo nodo del segmento s
                    r = self.nods.r[n]  # coordenadas de ese nodo
                    xx[f].append(r[0])
                    yy[f].append(r[1])
                grafs.append(
                    ax.plot(xx[f], yy[f], linestyle="-", marker="", label=str(f), color=sm.to_rgba(nc - 1 - c)))
        sm._A = []
        fig.colorbar(sm)

    @staticmethod
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    def pre_graficar_fibras(self, fig, ax, 
                            ncapas=None, 
                            lamr_min=None, lamr_max=None, 
                            color_por="nada",
                            byn=False, cmap='jet', barracolor=True,
                            colores_cm=None, ncolores_cm=20,
                            linewidth=2):
        lw = linewidth
        # preparo un mapa de colores mapeable por escalar
        if byn:
            mi_colormap = cm.get_cmap('gray_r')
            # lo trunco para que no incluya el cero (blanco puro que no hace contraste con el fondo)
            mi_colormap = self.truncate_colormap(mi_colormap, 0.2, 0.7)
        elif colores_cm is not None:
            mi_colormap = colors.LinearSegmentedColormap.from_list("mi_colormap", colores_cm, N=ncolores_cm)
        else:
            mi_colormap = cm.get_cmap(cmap)  # default: jet
            if colores_cm is not None:
                mi_colormap = colors.LinearSegmentedColormap.from_list("mi_colormap", colores_cm, N=ncolores_cm)
        if color_por == "lamr":
            # calculo los enrulamientos
            lamsr = self.calcular_enrulamientos()
            if lamr_min is None:
                lamr_min = np.min(lamsr)
            if lamr_max is None:
                lamr_max = np.max(lamsr)
            sm = plt.cm.ScalarMappable(cmap=mi_colormap, norm=plt.Normalize(vmin=lamr_min, vmax=lamr_max))
        elif color_por == "fibra":
            sm = plt.cm.ScalarMappable(cmap=mi_colormap, norm=plt.Normalize(vmin=0, vmax=len(self.fibs.con) - 1))
        elif color_por == "capa":
            sm = plt.cm.ScalarMappable(cmap=mi_colormap, norm=plt.Normalize(vmin=0, vmax=len(self.caps.con) - 1))
        elif color_por == "angulo":
            sm = plt.cm.ScalarMappable(cmap=mi_colormap, norm=plt.Normalize(vmin=0, vmax=PI))
        # dibujo las fibras (los segmentos)
        # preparo las listas, una lista para cada fibra
        xx = [list() for f in self.fibs.con]
        yy = [list() for f in self.fibs.con]
        grafs = list()
        if ncapas is None:
            caps_con = self.caps.con
        else:
            caps_con = self.caps.con[:ncapas]
        for c, c_con in enumerate(caps_con):  # recorro las capas
            for f in c_con:  # recorro las fibra de la capa
                f_con = self.fibs.con[f]
                # antes de recorrer los segmentos de cada fibra
                # el primer nodo del primer segmento lo agrego antes del bucle
                s = f_con[0]  # primer segmento de la fibra f
                n = self.segs.con[s][0]  # primer nodo del segmento s
                r = self.nods.r[n]  # coordenadas de ese nodo
                xx[f].append(r[0])
                yy[f].append(r[1])
                for s in f_con:  # recorro los segmentos de la fibra f
                    s_con = self.segs.con[s]
                    n = s_con[1]  # ultimo nodo del segmento s
                    r = self.nods.r[n]  # coordenadas de ese nodo
                    xx[f].append(r[0])
                    yy[f].append(r[1])
                if color_por == "lamr":
                    col = sm.to_rgba(lamsr[f])
                elif color_por == "fibra":
                    col = sm.to_rgba(f)
                elif color_por == "capa":
                    col = sm.to_rgba(c)
                    # aux = float(len(self.caps.con) - 1)
                    # lw = linewidth*(1. - 0.5*float(c)/aux)
                    # lw = linewidth*(0.7 + 0.3*float(c)/aux)
                elif color_por == "angulo":
                    theta = self.calcular_orientacion_extremo_extremo_de_una_fibra(f)
                    col = sm.to_rgba(theta)
                elif color_por == "nada":
                    col = "k"
                # print lw
                grafs.append(ax.plot(xx[f], yy[f], linestyle="-", marker="", label=str(f), color=col, linewidth=lw))
        if barracolor and color_por not in ("nada", "fibra"):
            sm._A = []
            cbar = fig.colorbar(sm)
            if color_por == "capa":
                cbar.set_ticks(range(len(self.caps.con)))

    def pre_graficar_interfibras(self, fig, ax, 
                                lamr_min=None, lamr_max=None, 
                                color_por="nada",
                                byn=False, colormap="jet", barracolor=True,
                                colores_cm=None, ncolores_cm=20):
        # preparo un mapa de colores mapeable por escalar
        # calculo las interfibras (conectividad)
        infbs_con = self.calcular_conectividad_de_interfibras()
        if byn:
            # mi_colormap = plt.cm.gray_r
            mi_colormap = cm.get_cmap('gray_r')
            # lo trunco para que no incluya el cero (blanco puro que no hace contraste con el fondo)
            mi_colormap = self.truncate_colormap(mi_colormap, 0.4, 1.0)
        elif colores_cm is not None:
                mi_colormap = colors.LinearSegmentedColormap.from_list("mi_colormap", colores_cm, N=ncolores_cm)
        else:
            mi_colormap = cm.get_cmap(colormap)
        if color_por == "lamr":
            # calculos los enrulamientos
            lamsr = self.calcular_enrulamientos_de_interfibras()
            if lamr_min is None:
                lamr_min = np.min(lamsr)
            if lamr_max is None:
                lamr_max = np.max(lamsr)
            sm = plt.cm.ScalarMappable(cmap=mi_colormap, norm=plt.Normalize(vmin=lamr_min, vmax=lamr_max))
        elif color_por == "interfibra":
            sm = plt.cm.ScalarMappable(cmap=mi_colormap, norm=plt.Normalize(vmin=0, vmax=len(infbs_con) - 1))
        elif color_por == "nada":
            sm = plt.cm.ScalarMappable(cmap=mi_colormap, norm=plt.Normalize(vmin=0, vmax=len(infbs_con) - 1))
        # dibujo las fibras (los segmentos)
        # preparo las listas, una lista para cada fibra
        xx = [list() for infb_con in infbs_con]
        yy = [list() for infb_con in infbs_con]
        grafs = list()
        print("-")
        print("graficando interfibras")
        nif = float(len(infbs_con))
        pctp = np.arange(0., 101., 10.).tolist()
        pcpd = np.zeros(len(pctp), dtype=bool).tolist()
        for i, infb_con in enumerate(infbs_con):  # recorro las interfibras
            # preparo algo para ir imprimiendo el porcentaje de progreso
            pc = round(float(i) / nif * 100., 0) 
            if pc in pctp:
                ipc = pctp.index(pc)
                if not pcpd[ipc]:
                    print("{:4.0f}% ".format(pc), end='')
                    pcpd[ipc] = True
            # antes de recorrer los segmentos de cada interfibra
            # el primer nodo del primer segmento lo agrego antes del bucle
            s = infb_con[0]  # primer segmento de la interfibra i
            n = self.segs.con[s][0]  # primer nodo del segmento s
            r = self.nods.r[n]  # coordenadas de ese nodo
            xx[i].append(r[0])
            yy[i].append(r[1])
            for s in infb_con:  # recorro los segmentos de la interfibra i
                s_con = self.segs.con[s]
                n = s_con[1]  # ultimo nodo del segmento s
                r = self.nods.r[n]  # coordenadas de ese nodo
                xx[i].append(r[0])
                yy[i].append(r[1])
            if color_por == "lamr":
                col = sm.to_rgba(lamsr[i])
            elif color_por == "interfibra":
                col = sm.to_rgba(i)
            elif color_por == "nada":
                col = "k"
            grafs.append(ax.plot(xx[i], yy[i], linestyle="-", marker="", label=str(i), color=col))
        print("-")
        if not byn and barracolor:
            sm._A = []
            fig.colorbar(sm)

    def pre_graficar_nodos_frontera(self, fig, ax, markersize=8):
        # dibujo las fibras (los segmentos)
        # preparo las listas, una lista para cada fibra
        xx = [list() for f in self.fibs.con]
        yy = [list() for f in self.fibs.con]
        # grafs = list() # un plot para cada fibra
        for f in range(len(self.fibs.con)):  # f es un indice
            # el primer nodo y el ultimo de cada fibra son fronteras
            s = self.fibs.con[f][0]  # obtengo el indice del primer segmento de la fibra numero f
            n = self.segs.con[s][0]  # obtengo el indice del primer nodo del segmento numero s
            r = self.nods.r[n]  # obtengo las coordenadas del nodo numero n
            xx[f].append(r[0])
            yy[f].append(r[1])
            s = self.fibs.con[f][-1]  # obtengo el indice del ultimo segmento de la fibra numero f
            n = self.segs.con[s][1]  # obtengo el indice del segundo nodo del ultimo numero s
            r = self.nods.r[n]  # obtengo las coordenadas del nodo numero n
            xx[f].append(r[0])
            yy[f].append(r[1])
            # grafs.append( ax.plot(xx[f], yy[f], linewidth=0, marker="x", mec="k", markersize=markersize) )
        ax.plot(xx, yy, linewidth=0, marker="o", mec="k", mfc="w", markersize=markersize)

    def pre_graficar_nodos_interseccion(self, fig, ax, markersize=8):
        # dibujo las fibras (los segmentos)
        # preparo las listas, una lista para cada fibra
        xx = []
        yy = []
        for n in range(len(self.nods.r)):
            if self.nods.tipos[n] == 2:
                xx.append(self.nods.r[n][0])
                yy.append(self.nods.r[n][1])
        ax.plot(xx, yy, linewidth=0, marker="o", mec="k", mfc="w", markersize=markersize)

    def pre_graficar_nodos_internos(self, fig, ax):
        # dibujo las fibras (los segmentos)
        # preparo las listas, una lista para cada fibra
        xx = []
        yy = []
        for n in range(len(self.nods.r)):
            if self.nods.tipos[n] == 0:
                xx.append(self.nods.r[n][0])
                yy.append(self.nods.r[n][1])
        ax.plot(xx, yy, linewidth=0, marker=".", markersize=1)

    def pre_graficar(self, fig, ax, lamr_min=None, lamr_max=None, byn=False):
        self.marco.graficar(fig, ax)
        self.pre_graficar_nodos_frontera(fig, ax)
        self.pre_graficar_nodos_interseccion(fig, ax)
        self.pre_graficar_nodos_internos(fig, ax)
        self.pre_graficar_fibras(fig, ax, lamr_min=lamr_min, lamr_max=lamr_max, byn=byn)
        # ax.legend(loc="upper left", numpoints=1, prop={"size":6})

    def graficar(self, fig=None, ax=None, lamr_min=None, lamr_max=None, byn=False):
        if ax is None:
            fig, ax = plt.subplots()
        self.pre_graficar(fig, ax, lamr_min, lamr_max, byn)
        plt.show()
