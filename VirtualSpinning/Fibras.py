class Fibras(object):
    """ Objeto de conjunto de fibras
     tiene 2 listas:
     :segs: conectividad como lista de listas de segmentos (solo los indices)
     :diams: lista de diametros con un valor por cada fibra
     """
    def __init__(self):
        self.segs = []  # conectividad (lista de listas de indices de segmentos)
        self.diams = []
        # self.l0s = []
        # self.th = []

    def add_fibra(self, segs, diam):
        self.segs.append(segs)
        self.diams.append(diam)

    def insertar_segmento(self, j, k, s):
        """ inserta un segmento en la conectividad de una fibra
        j: indice de la fibra
        k: indice donde se inserta el nuevo segmento
        s: indice del nuevo segmento para agregar a la conectividad """
        self.segs[j].insert(k, s)
