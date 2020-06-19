import numpy as np


SIDES = {
    'bottom': 0,
    'right': 1,
    'top': 2,
    'left': 3
}


class Marco(object):
    """ clase para armar el marco de la malla """ 
    def __init__(self, L):
        self.L = L
        self.nodes = np.array(
            [
                [0., 0.],
                [1., 0.],
                [1., 1.],
                [0., 1.]
            ]
        ) * L 
        self.sides = np.array(
            [
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0]
            ]
        )
        aux = self.nodes[self.sides]  # magia
        self.vecs = aux[:,1] - aux[:,0]
        self.uvecs = self.vecs / self.L

    def get_side_nodes(self, i):
        if isinstance(i, str):
            i = SIDES[i]
        return self.nodes[self.sides[i]] # ejemplo: i=1, self.nodes[[0,1]] -> [[0., 0.], [0., 1.]]

    def get_punto_random(self):
        boundary = np.random.randint(4)
        d = np.random.rand() * self.L
        r = self.nodes[boundary] + d * self.uvecs[boundary]
        return r, boundary     

    def graficar(self, fig, ax):
        # seteo limites
        margen = 0.1 * self.L 
        ax.set_xlim(left=0 - margen, right=self.L + margen)
        ax.set_ylim(bottom=0 - margen, top=self.L + margen)
        # dibujo los bordes del rve
        i = [0, 1, 2, 3, 0]
        x = self.nodes[i, 0]
        y = self.nodes[i, 1]
        ax.plot(x, y, linestyle=':', c='gray')