from VirtualSpinning.Aux.aux import append_to_keys


class Capas(object):
    """ Conjunto de capas donde cada capa es
     una conectividad de fibras """
    def __init__(self):
        self.num = 0
        self.con = []
        self.conT = {}  # conectividad traspuesta: fibras por capa

    def __getitem__(self, item):
        return self.con[item]

    def set_capas_listoflists(self, capas_con):
        self.num = len(capas_con)
        self.con = capas_con
        for j, capa_con in enumerate(capas_con):
            append_to_keys(self.conT, capa_con, j)

    def add_capa(self, cap_con):
        self.num += 1
        self.con.append(cap_con)
        append_to_keys(self.conT, cap_con, self.num - 1)

    def add_fibra_a_capa(self, jcapa, jfibra):
        """ jcapa y jfibra son los indices """
        self.con[jcapa].append[jfibra]