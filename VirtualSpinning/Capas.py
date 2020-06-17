class Capas(object):
    """ Conjunto de capas donde cada capa es
     una conectividad de fibras """
    def __init__(self):
        self.con = []

    def __getitem__(self, item):
        return self.con[item]

    def set_capas_listoflists(self, capas_con):
        self.con = capas_con

    def add_capa(self, cap_con):
        self.con.append(cap_con)

    def add_fibra_a_capa(self, jcapa, jfibra):
        """ jcapa y jfibra son los indices """
        self.con[jcapa].append[jfibra]