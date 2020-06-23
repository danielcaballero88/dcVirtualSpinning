from VirtualSpinning.aux import append_to_keys


class Fibras(object):
    """ Objeto de conjunto de fibras
     tiene 2 listas:
     :segs: conectividad como lista de listas de segmentos (solo los indices)
     :diams: lista de diametros con un valor por cada fibra
     """

    def __init__(self):
        self.num = 0
        self.con = []  # conectividad: va a ser una lista de listas de segmentos (sus indices nada mas)
        self.conT = {} # conectividad traspuesta
        self.dl = []  # longitud de segmento de cada fibra
        self.D = []  # diametro de cada fibra
        self.dth = []  # angulo de desviacion entre segmentos maximo de cada fibra
        self.loco = []  # longitud de contorno de la fibra

    def add_fibra(self, fib_con, dl, d, dtheta, loco):
        self.num += 1
        self.con.append(fib_con)
        append_to_keys(self.conT, fib_con, self.num-1)
        self.dl.append(dl)
        self.D.append(d)
        self.dth.append(dtheta)
        self.loco.append(loco)

    def actualizar_fibra(self, fib, loco):
        """
        Al mover un nodo es necesario actualizar algunos datos
        por lo pronto solamente la longitud de contorno (self.loco)
        """
        self.loco = loco