class Nodos(object):
    def __init__(self):
        self.num = 0
        self.r = []  # coordenadas de los nodos
        self.tipos = []  # lista de tipos (0=cont, 1=fron, 2=inter)

    def __getitem__(self, item):
        return self.r[item]

    def add_nodo(self, r_nodo, tipo):
        self.num += 1
        self.r.append(r_nodo)
        self.tipos.append(tipo)

    @staticmethod
    def append_to_key(dic, key, val):
        if key in dic: 
            dic[key].append(val) 
        else: 
            dic[key] = [val]

    def get_r(self, i_nodo):
        return self.r[i_nodo]

    def __len__(self):
        return len(self.r)