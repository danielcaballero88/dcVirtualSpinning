"""
Este modulo es para proveer funcionalidad de muestreo de
funciones estadistica. Dada una distribucion, se la inicia
con una larga lista de valores muestreados y luego al ir
necesitando valores, se los va sacando de esa lista
"""


import numpy as np
import scipy.stats as stats

# Numero por defecto para generar las distribuciones
BIGNUM = 100_000


class Rand(object):
    """
    Distribucion aleatoria entre dos valores
    """
    def __init__(self, n=BIGNUM, lower=0., upper=1.):
        self.n = n
        self.lower = lower
        self.upper = upper
        self.range = upper - lower
        self.x = self.lower + np.random.rand(n) * self.range
        self.i = 0

    def __call__(self):
        x = self.x[self.i]
        if self.i == self.n -1:
            self.i = 0
        else: 
            self.i += 1
        return x


class NormTr(object):
    def __init__(self, n=100000, loc=0., scale=1., lower=-10., upper=10.):
        self.n = n
        self.loc = loc
        self.scale = scale
        self.lower = (lower - loc) / scale
        self.upper = (upper - loc) / scale
        self.x = stats.truncnorm.rvs(self.lower, self.upper, loc=loc, scale=scale, size=self.n)
        self.i = 0

    def __call__(self):
        x = self.x[self.i]
        if self.i == self.n - 1:
            self.i = 0
        else:
            self.i += 1
        return x

    # def pdf(self, vec_x):
    #     return stats.truncnorm.rvs(self.lower, self.upper, loc=loc, scale=scale, size=self.n)


class Vals(object):
    def __init__(self, valores):
        self.n = len(valores)
        self.vals = valores
        self.i = 0

    def __call__(self):
        x = self.vals[self.i]
        if self.i == self.n - 1:
            self.i = 0
        else:
            self.i += 1
        return x
