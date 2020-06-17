class Nodos(object):
    def __init__(self):
        self.r = []  # coordenadas de los nodos
        self.tipos = []  # lista de tipos (0=cont, 1=fron, 2=inter)

    def __getitem__(self, item):
        return self.r[item]

    def add_nodo(self, r_nodo, tipo):
        self.r.append(r_nodo)
        self.tipos.append(tipo)

    def get_r(self, i_nodo):
        return self.r[i_nodo]

    def __len__(self):
        return len(self.r)