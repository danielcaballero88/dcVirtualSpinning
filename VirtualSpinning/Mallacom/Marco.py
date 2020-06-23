import numpy as np


SIDES = {
    'bottom': 0,
    'right': 1,
    'top': 2,
    'left': 3
}


MARCO_UNITARIO = np.array(
    [
        [0., 0.],
        [1., 0.],
        [1., 1.],
        [0., 1.]
    ]
)


class Marco(object):
    """ clase para armar el marco de la malla """ 
    def __init__(self, L):
        # Side length
        self.L = L
        # Coordinates
        self.nodes = MARCO_UNITARIO * L - 0.5*L
        self.left = self.nodes[:,0].min() 
        self.right = self.nodes[:,0].max() 
        self.top = self.nodes[:,1].max() 
        self.bottom = self.nodes[:,1].min()
        # Segments
        self.sides = np.array(
            [
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0]
            ]
        )
        # Vectors
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

    def check_punto_fuera(self, r):
        """
        metodo para chequear si r cae fuera del marco
        """
        x,y = r
        if x <= self.left or x >= self.right or y <= self.bottom or y >= self.top:
            return True
        else:
            return False

    def graficar(self, fig, ax):
        # seteo limites
        margen = 0.1 * self.L 
        ax.set_xlim(left=self.left - margen, right=self.right + margen)
        ax.set_ylim(bottom=self.bottom - margen, top=self.top + margen)
        # dibujo los bordes del rve
        i = [0, 1, 2, 3, 0]
        x = self.nodes[i, 0]
        y = self.nodes[i, 1]
        ax.plot(x, y, linestyle=':', c='gray')