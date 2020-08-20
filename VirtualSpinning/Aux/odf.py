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


def calcular_histograma(x, n=10, rango=None, polar=False):
    """ calcula el histograma de x dividiendo el intervalo
    en bins y contando el numero de ocurrencias en cada bin
    luego devuelve los centros de los bins, ancho (son todos iguales)
    y los valores de conteo
    para obtener frecuencias o porcentajes se hace por fuera de aqui"""
    x_count, x_edges = np.histogram(x, bins=n, range=rango)
    delta = (x_edges[-1] - x_edges[0]) / float(n)
    # pdf = conteo / float(np.sum(conteo)) / delta
    x_centers = x_edges[:-1] + 0.5*delta
    if polar:
        x_centers = np.concatenate((x_centers, x_centers[0:1]))
        x_count = np.concatenate((x_count, x_count[0:1]))
    return x_centers, delta, x_count


def calcular_histograma_relacionado(a, b, n=10, rango=None, polar=False):
    """ calcula la distribucion de b sobre un dominio dado por a
    se divide el intervalo de a en bins y se calcula, en cada bin
    el valor medio y la dispersion de los a"""
    # primero hago el histograma de a que es una gilada
    a_count, a_edges = np.histogram(a, bins=n, range=rango)
    a_delta = (a_edges[-1]-a_edges[0]) / float(n)
    a_centers = a_edges[:-1] + 0.5*a_delta
    # luego recorro los bins de a calculando dentro de cada uno
    # los valores o conteos de b que quiero
    b_means = np.zeros( (n,) , dtype=float )
    b_stds = np.zeros( (n,) , dtype=float )
    for j in range(n):
        # limites del bin
        a0,a1 = a_edges[j:j+2]
        # indices de los elementos de a que caen en el intervalo j
        jj = np.where( np.logical_and(a>=a0,a<a1) )
        # elementos de b que caen en el intervalo j de a
        bjj = b[jj]
        # ahora calculo el valor medio y lo guardo en el array
        b_means[j] = np.mean(bjj)
        # idem desv estandar
        b_stds[j] = np.std(bjj)
    # si la idea es hacerlo polar tengo que repetir al final el primer elemento
    if polar: 
        a_centers = np.concatenate((a_centers, a_centers[0:1]))
        b_means = np.concatenate((b_means, b_means[0:1]))
        b_stds = np.concatenate((b_stds, b_stds[0:1]))
    return a_centers, a_delta, b_means, b_stds


class DiscrDistr(object):
    """
    Base class for discrete distributions within this module 
    """ 
    def __init__(self, n, lower, upper): 
        """
        Args:
            n: number of points to sample
            lower: lower bound for the distribution 
            upper: upper bound for the distribution
        Computes:
            x: array of sampled values to withdraw numbers from
        """
        self.n = n 
        self.lower = lower 
        self.upper = upper
        self.range = upper - lower
        self.x = np.zeros(self.n)
        self.i = 0

    def __call__(self):
        """ 
        Each time the DiscrDistr is called, it return one value from its 
        array of sampled values, when hitting the end, it restarts
        """
        x = self.x[self.i]
        if self.i == self.n -1:
            self.i = 0
        else: 
            self.i += 1
        return x


class Rand(DiscrDistr):
    """
    Distribucion aleatoria entre dos valores
    """
    def __init__(self, n=BIGNUM, lower=0., upper=1.):
        super().__init__(n, lower, upper)
        self.x = self.lower + np.random.rand(n) * self.range


class NormTr(DiscrDistr):
    def __init__(self, n=BIGNUM, loc=0., scale=1., lower=-10., upper=10.):
        super().__init__(n, lower, upper)
        # Los valores que tengo que pasar al metodo de scipy no son exactamente lower y upper
        self.loc = loc
        self.scale = scale
        self.lowernt = (lower - loc) / scale
        self.uppernt = (upper - loc) / scale
        self.x = stats.truncnorm.rvs(self.lowernt, self.uppernt, loc=loc, scale=scale, size=self.n)


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
