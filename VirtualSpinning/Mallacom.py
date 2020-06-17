# Global Imports
import csv
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
# Local Imports
from VirtualSpinning.Capas import Capas
from VirtualSpinning.Fibras import Fibras
from VirtualSpinning.Segmentos import Segmentos
from VirtualSpinning.Nodos import Nodos
from VirtualSpinning.aux import iguales
from VirtualSpinning.aux import calcular_interseccion_entre_segmentos as calcular_interseccion
from VirtualSpinning.aux import find_string_in_file
from VirtualSpinning.aux import calcular_longitud_de_segmento
from VirtualSpinning.aux import calcular_angulo_de_segmento


class Mallacom(object):
    def __init__(self, L, Dm, volfrac, ls, devangmax, fundisor=None):
        self.L = L
        self.Dm = Dm  # diametro medio de fibras y espesor de las capas
        self.volfrac = volfrac  # si es float < 1 es volfrac, si es > 1 es num de fibras por capa
        self.ls = ls
        self.devangmax = devangmax
        self.fundisor = fundisor  # tiene que ser una funcion (or callable object) que devuelve un valor de orientacion
        self.caps = Capas()  # lista vacia
        self.fibs = Fibras()  # lista vacia
        self.segs = Segmentos()  # lista vacia
        self.nods = Nodos()  # tiene dos listas vacias
        self.bordes_n = Nodos()  # lista de coordenadas con los 4 nodos del borde
        self.bordes_s = Segmentos()  # lista con los segmentos con los 4 nodos del borde
        self.calcular_marco()
        self.pregraficado = False
        self.fig = None
        self.ax = None

    def calcular_marco(self):
        # agrego los 4 nodos
        self.bordes_n.add_nodo([0., 0.], 1)
        self.bordes_n.add_nodo([self.L, 0.], 1)
        self.bordes_n.add_nodo([self.L, self.L], 1)
        self.bordes_n.add_nodo([0., self.L], 1)
        # agrego los 4 segmentos
        self.bordes_s.add_segmento([0, 1], self.bordes_n.r)
        self.bordes_s.add_segmento([1, 2], self.bordes_n.r)
        self.bordes_s.add_segmento([2, 3], self.bordes_n.r)
        self.bordes_s.add_segmento([3, 0], self.bordes_n.r)

    def make_capa2(self, dl=None, d=None, dtheta=None, volfraction=None, orient_distr=None):
        """
        armo una capa con fibras, todas van a armarse con los
        mismos parmetros dl y dtheta (se debe modificar para usar distribuciones)
        se depositan fibras hasta que se supera la fraccion de volumen dictada
        """
        # si volfraction es int estoy dando el numero de fibras!
        if isinstance(volfraction, int):
            cond_fin_n = True
            n_final = volfraction
        else:
            cond_fin_n = False
            volc = self.L * self.L * self.Dm  # volumen de la capa
            vols_final = volfraction * volc  # volumen de solido (ocupado por fibras) a alcanzar
        # chequeo si uso los parametros globales de la malla o si di parametros diferentes para esta capa
        if dl is None:
            dl = self.ls
        if d is None:
            d = self.Dm
        if dtheta is None:
            dtheta = self.devangmax
        if volfraction is None:
            volfraction = self.volfrac
        if orient_distr is None:
            orient_distr = self.fundisor
        # --
        ncapas = len(self.caps.con)
        capa_con = list()
        i = 0
        vols = 0.  # volumen de solido actual
        while True:
            i += 1
            j = self.make_fibra2(dl, d, dtheta, orient_distr)
            if j == -1:
                i -= 1
            else:
                volf = self.calcular_volumen_de_una_fibra(j)
                vols += volf
                capa_con.append(j)
            # me fijo si complete la capa
            if cond_fin_n:
                if i == n_final:
                    break
            else:
                if vols >= vols_final:
                    break
        self.caps.add_capa(capa_con)

    def calcular_loco_de_una_fibra(self, f):
        """ calcula la longitud de contorno de una fibra """
        nsegs = len(self.fibs.con[f])
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
        dl = self.fibs.dls[f]
        d = self.fibs.ds[f]
        return loco * np.pi * d * d / 4.

    def calcular_fraccion_de_volumen_de_una_capa(self, capcon):
        """ calcula la fraccion de volumen de una capa como
        el volumen ocupado por fibras
        dividido el volumen total de la capa """
        volfs = 0
        for f in capcon:  # recorro las fibras de la capa
            volf = self.calcular_volumen_de_una_fibra(f)
            volfs += volf
        # el volumen total de la capa es:
        volc = self.L * self.L * self.Dm
        # luego la fraccion de volumen
        fracvol = volfs / volc
        return fracvol

    def calcular_orientacion_de_una_fibra(self, f):
        """ calcula la orientacion de una fibra como el promedio
        de las orientacions de sus segmentos
        cada orientacion es un angulo en [0,pi) """
        fcon = self.fibs.con[f]
        nsegs = len(fcon)
        theta_f = 0.
        for s in fcon:
            # s es un indice de segmento
            theta_s = self.segs.thetas[s]
            # theta_s esta en [0,pi)
            if theta_s >= np.pi:
                theta_s = theta_s - np.pi
            if theta_s < 0:
                pass
            # ahora voy haciendo el promedio
            theta_f += theta_s
        theta_f = theta_f / float(nsegs)
        return theta_f

    def calcular_orientacion_extremo_extremo_de_una_fibra(self, f):
        fcon = self.fibs.con[f]
        s0 = fcon[0]
        s1 = fcon[1]
        n0 = self.segs.con[s0][0]
        n1 = self.segs.con[s1][1]
        r0 = self.nods.r[n0]
        r1 = self.nods.r[n1]
        theta_2pi = calcular_angulo_de_segmento(r0, r1)
        if theta_2pi >= np.pi:
            theta = theta_2pi - np.pi
        else:
            theta = theta_2pi
        if theta_2pi < 0.:
            raise ValueError
        return theta

    def make_fibra2(self, dl, d, dtheta, orient_distr=None):
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
        coors = list()
        # Armo el primer segmento
        # primero busco un nodo en el contorno
        x0, y0, b0 = self.get_punto_sobre_frontera()
        if orient_distr is None:
            theta_abs = np.random.rand() * np.pi
        else:
            theta_abs = orient_distr()
        # theta_abs = 179. * np.pi/180.
        if theta_abs == np.pi:
            theta_abs = 0.
        elif theta_abs > np.pi:
            raise ValueError("theta_abs de una fibra no comprendido en [0,pi)")
        # veo el cuadrante
        if theta_abs < np.pi * 1.0e-8:
            cuad = -1  # direccion horizontal
        elif np.abs(theta_abs - np.pi * 0.5) < 1.0e-8:
            cuad = -2  # direccion vertical
        elif theta_abs < np.pi * 0.5:
            cuad = 1  # primer cuadrante
        else:
            cuad = 2  # segundo cuadrante
        # ahora me fijo la relacion entre cuadrante y borde
        if cuad == -1:  # fibra horizontal
            if b0 in (0, 2):
                return -1  # esta fibra no vale, es horizontal sobre un borde horizontal
            elif b0 == 1:
                theta = np.pi
            else:  # b0 == 3
                theta = 0.
        elif cuad == -2:  # fibra vertical
            if b0 in (1, 3):
                return -1
            elif b0 == 0:
                theta = 0.5 * np.pi
            else:  # b0 == 2
                theta = 1.5 * np.pi
        elif cuad == 1:  # primer cuadrante
            if b0 in (0, 3):
                theta = theta_abs
            else:  # b0 in(1,2)
                theta = theta_abs + np.pi
        else:  # cuad == 2 segundo cuadrante
            if b0 in (0, 1):
                theta = theta_abs
            else:  # b0 in (2,3)
                theta = theta_abs + np.pi
        # ya tengo el angulo del segmento
        dx = dl * np.cos(theta)
        dy = dl * np.sin(theta)
        coors.append([x0, y0])
        coors.append([x0 + dx, y0 + dy])
        # ahora agrego nuevos nodos en un bucle
        # cada iteracion corresponde a depositar un nuevo segmento
        n = 1
        while True:
            # si el nodo anterior ha caido fuera del rve ya esta la fibra
            if self.check_fuera_del_RVE(coors[-1]):
                break
            n += 1
            # de lo contrario armo un nuevo segmento a partir del ultimo nodo
            # el angulo puede sufrir variacion
            theta = theta + dtheta * (2.0 * np.random.rand() - 1.0)
            # desplazamiento:
            dx = dl * np.cos(theta)
            dy = dl * np.sin(theta)
            # nuevo nodo
            x = coors[-1][0] + dx
            y = coors[-1][1] + dy
            coors.append([x, y])
        # -
        # Aqui termine de obtener las coordenadas de los nodos que componen la fibra
        # si la fibra es muy corta la voy a descartar
        # para eso calculo su longitud de contorno
        loco = dl * float(len(coors) - 1)  # esto es aproximado porque el ultimo segmento se recorta
        if loco < 0.3 * self.L:
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
        self.trim_fibra_at_frontera(f_con)  # lo comento porque a veces quedan segmentos super pequenos
        self.fibs.add_fibra(f_con, dl, d, dtheta)
        return len(self.fibs.con) - 1  # devuelvo el indice de la fibra

    def get_punto_sobre_frontera(self):
        boundary = np.random.randint(4)
        d = np.random.rand() * self.L
        if boundary == 0:
            x = d
            y = 0.0
        elif boundary == 1:
            x = self.L
            y = d
        elif boundary == 2:
            x = self.L - d
            y = self.L
        elif boundary == 3:
            x = 0.0
            y = self.L - d
        return x, y, float(boundary)

    def check_fuera_del_RVE(self, r):
        x = r[0]
        y = r[1]
        if x <= 0 or x >= self.L or y <= 0 or y >= self.L:
            return True
        else:
            return False

    def trim_fibra_at_frontera(self, fib_con):
        """ subrutina para cortar la fibra que ha salido del rve """
        # debo cortar la ultima fibra en su interseccion por el rve
        # para eso calculo las intersecciones de los nodos con los bordes
        # coordenadas del ultimo segmento de la fibra de conectividad fib_con
        s = fib_con[-1]
        rs0 = self.nods.r[self.segs.con[s][0]]  # coordenadas xy del nodo 0 del segmento s
        rs1 = self.nods.r[self.segs.con[s][1]]  # coordenadas xy del nodo 1 del segmento s
        # pruebo con cada borde
        for b in range(len(self.bordes_s.con)):  # recorro los 4 bordes
            # puntos del borde en cuestion
            rb0 = self.bordes_n.r[self.bordes_s.con[b][0]]  # coordenadas xy del nodo 0 del borde b
            rb1 = self.bordes_n.r[self.bordes_s.con[b][1]]  # coordenadas xy del nodo 1 del borde b
            interseccion = calcular_interseccion(rs0, rs1, rb0, rb1)
            if interseccion is None:  # no hubo interseccion
                continue  # con este borde no hay interseccion, paso al que sigue
            else:  # hubo interseccion
                in_r, in_tipo, in_seg0, in_seg1 = interseccion
                if in_tipo == 2:  # interseccion en el medio
                    try:  # tengo que mover el ultimo nodo y por lo tanto cambia el segmento
                        self.segs.mover_nodo(s, 1, self.nods.r, in_r)
                        rs1 = in_r
                    except ValueError as e:
                        print("error")
                        print(fib_con, b, interseccion)
                        quit()
                else:  # interseccion coincide con uno o dos extremos
                    # en este caso solo me importa el segundo nodo del primer segmento (seg 0)
                    # porque el segmento 1 es el borde, y el primer nodo del seg 0 siempre deberia estar dentro del rve
                    # (o en el borde a lo sumo si se trata de una fibra de un solo segmento)
                    # y en ese caso no hay nada que hacer! puesto que el nodo ya esta en el borde
                    pass

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
        """ ojo son diferentes a las subfibras de una malla simplificada
        aqui las interfibras son concatenaciones de segmentos entre nodos interseccion
        en una ms las subfibras son una simplificacion de una fibra dando solamente los nodos extremos y enrulamiento """
        infbs_con = list()  # conectividad: lista de listas de segmentos
        for f, fcon in enumerate(self.fibs.con):  # recorro las fibras
            # cada fibra que empieza implica una nueva interfibra
            infb = list()  # conectividad de la interfibra: lista de segmentos
            # tengo que ir agregando segmentos hasta toparme con un nodo interseccion o frontera
            for s in fcon:  # recorro los segmentos de la fibra f
                scon = self.segs.con[s]
                n0, n1 = scon
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
        for f, fcon in enumerate(self.fibs.con):
            # theta_f = self.calcular_orientacion_de_una_fibra(f)
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
        conteo, x_edges = np.histogram(phis, bins=n, range=(0., np.pi))
        delta = (np.pi - 0.) / float(n)
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
        if lamr_min is None:
            lamr_min = np.min(lamsr)
        if lamr_max is None:
            lamr_max = np.max(lamsr)
        # me fijo si hay imposicion de ancho de bin, entonces calculo con eso el numero de bins
        if binwidth is not None:
            n = int((lamr_max - lamr_min) / binwidth + 0.5)  # es el entero mas cercano
        conteo, x_edges = np.histogram(lamsr, bins=n, range=(lamr_min, lamr_max))
        delta = (lamr_max - lamr_min) / n
        # pdf = conteo / float(np.sum(conteo)) / delta
        x = x_edges[:-1] + 0.5 * delta
        return x, delta, conteo

    def get_histograma_lamr(self, lamr_min=None, lamr_max=None, nbins=5, binwidth=None, opcion="fibras",
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

    def calcular_enrulamientos_de_interfibras(self):
        """ calcular para todas las interfibras sus longitudes de contorno y
        sus longitudes extremo a extremos (loco y lete)
        y calcula el enrulamiento como lamr=loco/lete """
        lamsr = []
        infbs_con = self.calcular_conectividad_de_interfibras()
        for infb_con in infbs_con:  # recorro las interfibras (fibras interectadas) del rve
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
        # primero escribo L, dl y dtheta
        fid.write("*Parametros (L, Dm, volfrac, ls, devangmax) \n")
        fid.write("{:20.8f}\n".format(self.L))
        fid.write("{:20.8f}\n".format(self.Dm))
        fid.write("{:20.8f}\n".format(self.volfrac))
        fid.write("{:20.8f}\n".format(self.ls))
        fid.write("{:20.8f}\n".format(self.devangmax * 180. / np.pi))  # lo escribo en grados
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
        # sigo con las fibras: indice, dl, d, dtheta, nsegs_f, y segmentos (conectividad)
        dString = "*Fibras \n" + str(len(self.fibs.con)) + "\n"
        fid.write(dString)
        for f, fcon in enumerate(self.fibs.con):
            dString = "{:12d}".format(f)  # indice
            dString += "{:17.8e}{:17.8e}{:17.8e}".format(self.fibs.dls[f], self.fibs.ds[f],
                                                         self.fibs.dthetas[f])  # dl, d y dtheta
            dString += "{:12d}".format(len(fcon))  # indice
            dString += "".join("{:12d}".format(val) for val in fcon) + "\n"  # conectividad
            fid.write(dString)
        # termino con las capas: indice y fibras (conectividad):
        dString = "*Capas \n" + str(len(self.caps.con)) + "\n"
        fid.write(dString)
        for c, ccon in enumerate(self.caps.con):
            dString = "{:12d}".format(c)  # indice
            dString += "{:12d}".format(len(ccon))  # indice
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
        ierr = find_string_in_file(fid, target, True)
        L = float(fid.next())
        Dm = float(fid.next())
        volfrac = float(fid.next())
        if volfrac > 1.:
            volfrac = int(volfrac + .5)
        ls = float(fid.next())
        devangmax = float(fid.next())  # en grados
        devangmax = devangmax * np.pi / 180.
        # luego busco coordenadas
        target = "*coordenadas"
        ierr = find_string_in_file(fid, target, True)
        num_r = int(fid.next())
        coors = list()
        tipos = list()
        for i in range(num_r):
            j, t, x, y = (float(val) for val in fid.next().split())
            tipos.append(int(t))
            coors.append([x, y])
        # luego los segmentos
        target = "*segmentos"
        ierr = find_string_in_file(fid, target, True)
        num_s = int(fid.next())
        segs = list()
        for i in range(num_s):
            j, n0, n1 = (int(val) for val in fid.next().split())
            segs.append([n0, n1])
        # luego las fibras
        target = "*fibras"
        ierr = find_string_in_file(fid, target, True)
        num_f = int(fid.next())
        fibs = list()
        dls = list()
        ds = list()
        dthetas = list()
        for i in range(num_f):
            svals = fid.next().split()
            j = int(svals[0])
            dl = float(svals[1])
            d = float(svals[2])
            dtheta = float(svals[3])
            nsegsf = int(svals[4])
            fcon = [int(val) for val in svals[5:]]
            fibs.append(fcon)
            dls.append(dl)
            ds.append(d)
            dthetas.append(dtheta)
        # luego la capas
        target = "*capas"
        ierr = find_string_in_file(fid, target, True)
        num_c = int(fid.next())
        caps = list()
        for c in range(num_c):
            svals = fid.next().split()
            j = int(svals[0])
            nfibsc = int(svals[1])
            ccon = [int(val) for val in svals[2:]]
            caps.append(ccon)
        # ahora que tengo todo armo el objeto
        malla = cls(L, Dm, volfrac, ls, devangmax)
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
            dl = dls[i]
            d = ds[i]
            dtheta = dthetas[i]
            malla.fibs.add_fibra(f_con, dl, d, dtheta)
        # le asigno las capas
        for c in range(num_c):
            c_con = caps[c]
            malla.caps.add_capa(c_con)
        # listo
        return malla

    def pre_graficar_bordes(self, fig, ax, byn=False):
        # seteo
        margen = 0.1 * self.L
        ax.set_xlim(left=0 - margen, right=self.L + margen)
        ax.set_ylim(bottom=0 - margen, top=self.L + margen)
        # dibujo los bordes del rve
        fron = []
        fron.append([[0, self.L], [0, 0]])
        fron.append([[0, 0], [self.L, 0]])
        fron.append([[0, self.L], [self.L, self.L]])
        fron.append([[self.L, self.L], [self.L, 0]])
        plt_fron0 = ax.plot(fron[0][0], fron[0][1], linestyle=":", c="gray")
        plt_fron1 = ax.plot(fron[1][0], fron[1][1], linestyle=":", c="gray")
        plt_fron2 = ax.plot(fron[2][0], fron[2][1], linestyle=":", c="gray")
        plt_fron3 = ax.plot(fron[3][0], fron[3][1], linestyle=":", c="gray")

    def pre_graficar_capas(self, fig, ax, byn=True):
        nc = len(self.caps.con)
        if byn:
            mi_colormap = plt.cm.gray
        else:
            mi_colormap = plt.cm.rainbow
        sm = plt.cm.ScalarMappable(cmap=mi_colormap, norm=plt.Normalize(vmin=0, vmax=nc - 1))
        # dibujo las fibras (los segmentos)
        # preparo las listas, una lista para cada fibra
        xx = [list() for f in self.fibs.con]
        yy = [list() for f in self.fibs.con]
        grafs = list()
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

    def pre_graficar_fibras(self, fig, ax, ncapas=None, lamr_min=None, lamr_max=None, byn=False, barracolor=True,
                            color_por="nada", linewidth=2, colores_cm=None, ncolores_cm=20):
        lw = linewidth
        # preparo un mapa de colores mapeable por escalar
        lamsr = self.calcular_enrulamientos()
        if byn:
            mi_colormap = plt.cm.gray_r
            # lo trunco para que no incluya el cero (blanco puro que no hace contraste con el fondo)
            mi_colormap = self.truncate_colormap(mi_colormap, 0.2, 0.7)
        else:
            mi_colormap = plt.cm.jet
            if colores_cm is not None:
                mi_colormap = colors.LinearSegmentedColormap.from_list("mi_colormap", colores_cm, N=ncolores_cm)
        if color_por == "lamr":
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
            sm = plt.cm.ScalarMappable(cmap=mi_colormap, norm=plt.Normalize(vmin=0, vmax=np.pi))
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

    def pre_graficar_interfibras(self, fig, ax, lamr_min=None, lamr_max=None, byn=False, barracolor=True,
                                 color_por="nada", colormap="jet", colores_cm=None, ncolores_cm=20):
        # preparo un mapa de colores mapeable por escalar
        infbs_con = self.calcular_conectividad_de_interfibras()
        lamsr = self.calcular_enrulamientos_de_interfibras()
        if byn:
            mi_colormap = plt.cm.gray_r
            # lo trunco para que no incluya el cero (blanco puro que no hace contraste con el fondo)
            mi_colormap = self.truncate_colormap(mi_colormap, 0.4, 1.0)
        else:
            if colores_cm is not None:
                mi_colormap = colors.LinearSegmentedColormap.from_list("mi_colormap", colores_cm, N=ncolores_cm)
            elif colormap == "jet":
                mi_colormap = plt.cm.jet
            elif colormap == "prism":
                mi_colormap = plt.cm.prism
            elif colormap == "Dark2":
                mi_colormap = plt.cm.Dark2
        if color_por == "lamr":
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
        pcc = 0.
        nif = float(len(infbs_con))
        pctp = np.arange(0., 101., 10.).tolist()
        pcpd = np.zeros(len(pctp), dtype=bool).tolist()
        for i, infb_con in enumerate(infbs_con):  # recorro las interfibras
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
        xx = list()
        yy = list()
        # grafs = list() # un plot para cada fibra
        for n in range(len(self.nods.r)):
            if self.nods.tipos[n] == 2:
                xx.append(self.nods.r[n][0])
                yy.append(self.nods.r[n][1])
        ax.plot(xx, yy, linewidth=0, marker="o", mec="k", mfc="w", markersize=markersize)

    def pre_graficar_nodos_internos(self, fig, ax):
        # dibujo las fibras (los segmentos)
        # preparo las listas, una lista para cada fibra
        xx = list()
        yy = list()
        grafs = list()  # un plot para cada fibra
        for n in range(len(self.nods.r)):
            if self.nods.tipos[n] == 0:
                xx.append(self.nods.r[n][0])
                yy.append(self.nods.r[n][1])
        ax.plot(xx, yy, linewidth=0, marker=".", markersize=1)

    def pre_graficar(self, fig, ax, lamr_min=None, lamr_max=None, byn=False):
        self.pre_graficar_bordes(fig, ax, byn)
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
