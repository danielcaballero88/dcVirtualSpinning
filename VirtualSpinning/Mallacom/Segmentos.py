import numpy as np
from VirtualSpinning.aux import iguales
from VirtualSpinning.aux import calcular_angulo_de_segmento
from VirtualSpinning.aux import append_to_keys

class Segmentos(object):
    """
    La clase Segmentos funciona como un conjunto de segmentos
    Tiene tres listas: coordenadas (r), angulos (th), longitudes iniciales (l0)
    """

    def __init__(self):
        self.num = 0
        self.con = []  # lista de listas de dos nodos (indices)
        self.conT = {} # dict de listas con los segmentos que tienen cada nodo
        self.thetas = []
        self.longs = []

    def __len__(self):
        return len(self.con)

    def add_segmento(self, seg_con, coors):
        """
        aca las coordenadas las necesito para calcularle a cada segmento su longitud y angulo
        seg_con es la conectividad (2 nodos) del segmento
        coors son las coordenadas (lista de listas de a dos floats) de todos los nodos
        (con todos los nodos hasta el momento de crear este segmento esta bien,
        alcanza con que esten presentes en la lista los dos nodos de seg_con)
        intersec indica si el segmento ha sido intersectado aun o no"""
        self.num += 1
        self.con.append(seg_con)
        append_to_keys(dic=self.conT, keys=seg_con, val=self.num-1)
        self.calc_long(j=-1, coors=coors, new=True)
        self.calc_theta(j=-1, coors=coors, new=True)

    def actualizar_segmento(self, j, coors):
        """ en caso de que se mueva un nodo y haya que actualizar theta y longitud """
        self.calc_long(j, coors) 
        self.calc_theta(j, coors)

    # def cambiar_conectividad(self, j, new_con, coors):
    #     """ se modifica la conectividad de un segmento (j) de la lista
    #     se le da la nueva conectividad new_con
    #     y por lo tanto se vuelve a calcular su angulo y longitud
    #     (util para dividir segmentos en 2) """
    #     # cambio la conectividad 
    #     self.con[j] = new_con
    #     # TODO: cambiar la conectividad traspuesta
    #     self.calc_long(j, coors) 
    #     self.calc_theta(j, coors)

    def calc_long(self, j, coors, new=False):
        """
        Calcular la longitud de un segmento, si es un segmento nuevo
        entonces se debe anexar a la lista de longitudes, 
        si es un segmento viejo se debe modificar su valor en la lista

        Args:
            j: indice del segmento
            coors: array de coordenadas de los nodos 
            new: boolean, True si es un segmento nuevo
        """
        n0, n1 = self.con[j]
        dr = coors[n1] - coors[n0]
        long = np.sqrt(np.sum(dr*dr))
        if new:
            self.longs.append(long)
        else: 
            self.longs[j] = long

    def calc_theta(self, j, coors, new=False):
        """
        Calcular el angulo de un segmento, si es un segmento nuevo
        entonces se debe anexar a la lista de angulos, 
        si es un segmento viejo se debe modificar su valor en la lista

        Args:
            j: indice del segmento
            coors: array de coordenadas de los nodos 
            new: boolean, True si es un segmento nuevo
        """
        n0, n1 = self.con[j] 
        r0, r1 = coors[[n0, n1]]
        theta = calcular_angulo_de_segmento(r0, r1)
        if new:
            self.thetas.append(theta)
        else: 
            self.thetas[j] = theta

    def get_right(self, j, coors):
        n0 = self.con[j][0]
        n1 = self.con[j][1]
        x0 = coors[n0][0]
        x1 = coors[n1][0]
        return np.maximum(x0, x1)

    def get_left(self, j, coors):
        n0 = self.con[j][0]
        n1 = self.con[j][1]
        x0 = coors[n0][0]
        x1 = coors[n1][0]
        return np.minimum(x0, x1)

    def get_top(self, j, coors):
        n0 = self.con[j][0]
        n1 = self.con[j][1]
        y0 = coors[n0][1]
        y1 = coors[n1][1]
        return np.maximum(y0, y1)

    def get_bottom(self, j, coors):
        n0 = self.con[j][0]
        n1 = self.con[j][1]
        y0 = coors[n0][1]
        y1 = coors[n1][1]
        return np.minimum(y0, y1)

    def get_dx(self, j, coors):
        n0 = self.con[j][0]
        n1 = self.con[j][1]
        x0 = coors[n0][0]
        x1 = coors[n1][0]
        return x1 - x0

    def get_dy(self, j, coors):
        n0 = self.con[j][0]
        n1 = self.con[j][1]
        y0 = coors[n0][1]
        y1 = coors[n1][1]
        return y1 - y0

    def get_dx_dy_brtl(self, j, coors):
        n0 = self.con[j][0]
        n1 = self.con[j][1]
        x0 = coors[n0][0]
        y0 = coors[n0][1]
        x1 = coors[n1][0]
        y1 = coors[n1][1]
        return x1 - x0, y1 - y0, np.minimum(y0, y1), np.maximum(x0, x1), np.maximum(y0, y1), np.minimum(x0, x1)
