class Fibras(object):
    """ Objeto de conjunto de fibras
     tiene 2 listas:
     :segs: conectividad como lista de listas de segmentos (solo los indices)
     :diams: lista de diametros con un valor por cada fibra
     """

    def __init__(self):
        self.con = []  # conectividad: va a ser una lista de listas de segmentos (sus indices nada mas)
        self.dls = []
        self.ds = []
        self.dthetas = []

    def add_fibra(self, fib_con, dl, d, dtheta):
        self.con.append(fib_con)
        self.dls.append(dl)
        self.ds.append(d)
        self.dthetas.append(dtheta)

    def insertar_segmento(self, j, k, s):
        """ inserta un segmento en la conectividad de una fibra
        j: indice de la fibra
        k: indice donde se inserta el nuevo segmento
        s: indice del nuevo segmento para agregar a la conectividad """
        self.con[j].insert(k, s)
