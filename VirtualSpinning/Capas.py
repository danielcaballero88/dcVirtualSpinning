class Capas(object):
    """ Conjunto de capas donde cada capa es
     una conectividad de fibras """
    def __init__(self):
        self.fibs = []

    def __getitem__(self, item):
        return self.fibs[item]

    def set_capas_listoflists(self, capas_con):
        self.fibs = capas_con

    def add_capa(self, cap_con):
        self.fibs.append(cap_con)

    def add_fibra_a_capa(self, jcapa, kfibra):
        """ jcapa y kfibra son los indices """
        self.fibs[jcapa].append[kfibra]