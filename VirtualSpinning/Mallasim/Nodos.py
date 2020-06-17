import numpy as np


class Nodos(object):
    def __init__(self, n, r0, tipos):
        self.n = n
        self.r0 = np.array(r0, dtype=float)
        self.r = self.r0.copy()
        self.t = np.array(tipos, dtype=int)
        self.mf = self.t == 1
        self.mi = self.t == 2
