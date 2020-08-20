import numpy as np


SIDES = {
    'bottom': 0,
    'right': 1,
    'top': 2,
    'left': 3
}


MARCO_UNITARIO = np.array(
    [
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0]
    ]
)


IDMAT = np.array(
    [
        [1.0, 0.0],
        [0.0, 1.0]
    ]
)


class Marco(object):
    """ clase para armar el marco de la malla """ 
    def __init__(self, L):
        # Side length
        self.L = L
        # Coordinates
        self.nodes = MARCO_UNITARIO * L - 0.5*L
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

    def get_bottom(self):
        return self.nodes[:,1].min() 

    def get_right(self):
        return self.nodes[:,0].max() 

    def get_top(self):
        return self.nodes[:,1].max() 

    def get_left(self):
        return self.nodes[:,0].min() 
        

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
        outwards = (
            y <= self.get_bottom(),
            y >= self.get_top(),
            y >= self.get_top(),
            x <= self.get_left()
        )
        if any(outwards):
            return True
        else:
            return False

    def check_puntos_fuera(self, r):
        """
        Args:
            r: np.ndarray de posiciones
        Result:
            mask_out: tuple de mask arrays indicando True donde el nodo cae
                      fuera del rve en el orden: bottom, right, top, left
        """
        x = r[:,0]
        y = r[:,1]
        mask_bottom = y < self.get_bottom() 
        mask_right = x > self.get_right()
        mask_top = y > self.get_top()
        mask_left = x < self.get_left()
        return mask_bottom, mask_right, mask_top, mask_left

    def deformar(self, Fmacro):
        # coordenadas deformadas
        self.nodes = np.matmul(self.nodes, np.transpose(Fmacro))

    def graficar(self, fig, ax, limites={}):
        # seteo limites
        margen = 0.1 * self.L 
        bottom = limites.get('bottom', self.get_bottom() - margen)
        right = limites.get('right', self.get_right() + margen)
        top = limites.get('top', self.get_top() + margen)
        left = limites.get('left', self.get_left() - margen)
        ax.set_xlim(left=left, right=right)
        ax.set_ylim(bottom=bottom, top=top)
        # dibujo los bordes del rve
        i = [0, 1, 2, 3, 0]
        x = self.nodes[i, 0]
        y = self.nodes[i, 1]
        ax.plot(x, y, linestyle=':', c='gray')