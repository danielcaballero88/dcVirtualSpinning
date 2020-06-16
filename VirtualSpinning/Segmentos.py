import numpy as np
from VirtualSpinning.aux import iguales


class Segmentos(object):
    """
    La clase Segmentos funciona como un conjunto de segmentos
    Tiene tres listas: coordenadas (r), angulos (th), longitudes iniciales (l0)
    """

    def __init__(self):
        self.num = 0
        self.nods = []  # nodos (conectividad)
        self.ths = []  # angulos
        self.l0s = []  # longitudes iniciales

    def __len__(self):
        """
        para poder hacer len(self)
        """
        return self.num

    def add_segmento(self, nods, coords):
        """
        agrego un segmento al conjunto
        """
        self.nods.append(nods)
        try:
            longitud, angulo = self.calcular_long_y_theta(nods, coords)
        except ValueError:
            raise ValueError("Error, segmento de longitud nula!!")
        self.ths.append(angulo)
        self.l0s.append(longitud)

    def actualizar_segmento(self, j, coors):
        """ en caso de que se mueva un nodo y haya que actualizar theta y longitud """
        long, ang = self.calcular_long_y_theta(self.nods[j], coors)
        self.ths[j] = ang
        self.l0s[j] = long

    def mover_nodo(self, j, n, coors, new_r):
        """ mueve un nodo del segmento
        coors es una lista, es un objeto mutable
        por lo que al salir de este metodo se va ver modificada
        es decir, es un puntero
        j es el indice del segmento a moverle un nodo
        n es el indice del nodo para el segmento: 0 es inicial, 1 es final """
        assert n in (0, 1)
        nglobal = self.nods[j][n]
        coors[nglobal] = new_r  # se lo modifica resida donde resida (normalmente en un objeto nodos)
        self.actualizar_segmento(j, coors)

    def cambiar_conectividad(self, j, new_con, coors):
        """ se modifica la conectividad de un segmento (j) de la lista
        se le da la nueva conectividad new_con
        y por lo tanto se vuelve a calcular su angulo y longitud
        (util para dividir segmentos en 2) """
        self.nods[j] = new_con
        longitud, angulo = self.calcular_long_y_theta(new_con, coors)
        self.ths[j] = angulo
        self.l0s[j] = longitud

    @staticmethod
    def calcular_long_y_theta(seg, coors):
        n0 = seg[0]
        n1 = seg[1]
        dx = coors[n1][0] - coors[n0][0]
        dy = coors[n1][1] - coors[n0][1]
        long = np.sqrt(dx * dx + dy * dy)
        # ahora theta
        if iguales(dx, 0.0):
            # segmento vertical
            if iguales(dy, 0.0, 1.0e-12):
                raise ValueError("Error, segmento de longitud nula!!")
            elif dy > 0:
                theta = np.pi * .5
            else:
                theta = 1.5 * np.pi
        elif iguales(dy, 0):
            # segmento horizontal
            if dx > 0:
                theta = 0.0
            else:
                theta = np.pi
        else:
            # segmento oblicuo
            if dx < 0:
                # segundo o tercer cuadrante
                theta = np.pi + np.arctan(dy / dx)
            elif dy > 0:
                # primer cuadrante (dx>0)
                theta = np.arctan(dy / dx)
            else:
                # dx>0 and dy<0
                # cuarto cuadrante
                theta = 2.0 * np.pi + np.arctan(dy / dx)
        return long, theta

    def get_right(self, j, coors):
        n0 = self.nods[j][0]
        n1 = self.nods[j][1]
        x0 = coors[n0][0]
        x1 = coors[n1][0]
        return np.maximum(x0, x1)

    def get_left(self, j, coors):
        n0 = self.nods[j][0]
        n1 = self.nods[j][1]
        x0 = coors[n0][0]
        x1 = coors[n1][0]
        return np.minimum(x0, x1)

    def get_top(self, j, coors):
        n0 = self.nods[j][0]
        n1 = self.nods[j][1]
        y0 = coors[n0][1]
        y1 = coors[n1][1]
        return np.maximum(y0, y1)

    def get_bottom(self, j, coors):
        n0 = self.nods[j][0]
        n1 = self.nods[j][1]
        y0 = coors[n0][1]
        y1 = coors[n1][1]
        return np.minimum(y0, y1)

    def get_dx(self, j, coors):
        n0 = self.nods[j][0]
        n1 = self.nods[j][1]
        x0 = coors[n0][0]
        x1 = coors[n1][0]
        return x1 - x0

    def get_dy(self, j, coors):
        n0 = self.nods[j][0]
        n1 = self.nods[j][1]
        y0 = coors[n0][1]
        y1 = coors[n1][1]
        return y1 - y0

    def get_dx_dy_brtl(self, j, coors):
        n0 = self.nods[j][0]
        n1 = self.nods[j][1]
        x0 = coors[n0][0]
        y0 = coors[n0][1]
        x1 = coors[n1][0]
        y1 = coors[n1][1]
        return x1 - x0, y1 - y0, np.minimum(y0, y1), np.maximum(x0, x1), np.maximum(y0, y1), np.minimum(x0, x1)
