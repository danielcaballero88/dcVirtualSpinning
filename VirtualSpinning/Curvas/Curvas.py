import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import warnings

SMALL_SIZE = 12
MEDIUM_SIZE = 20
BIGGER_SIZE = 24
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

class Graficador(object):
    def __init__(self):
        self.d = {}

class Curva(object):
    def __init__(self, nombre, x, y):
        self.nombre = nombre
        self.n = x.shape[0]
        assert y.shape[0] == self.n
        self.x = x
        self.y = y
        self.d = None # derivada

    @classmethod
    def leer_de_csv(cls, nombre, csv_file, xcol=0, ycol=1, **kwargs):
        df = pd.read_csv(csv_file, **kwargs)
        x = df.iloc[:,xcol].values
        y = df.iloc[:,ycol].values
        c = cls(nombre, x, y)
        return c

    def get_copia(self, nombre):
        c = self.__class__(nombre, self.x.copy(), self.y.copy())
        return c

    def reducir_segun_dx(self, dx):
        x = np.arange(self.x[0],self.x[-1],dx)
        y = np.interp(x, self.x, self.y)
        self.x = x
        self.y = y

    def reducir_segun_dx_promediando(self, dx):
        x = np.arange(self.x[0],self.x[-1]-.5*dx,dx) + .5*dx
        y = np.zeros(shape=x.shape, dtype=float)
        for i,xi in enumerate(x):
            yi = self.y[np.logical_and(self.x>xi-.5*dx, self.x<xi+.5*dx)]
            if yi.size == 0:
                warnings.warn('dx muy pequeno!')
                return
            else:
                y[i] = yi.mean()
        self.x = x
        self.y = y

    def reducir_segun_x_promediando(self, x):
        y = np.zeros(shape=x.shape, dtype=float)
        dx = x[1:] - x[:-1] # ojo que tiene dimension menor a x en 1 item
        dx = np.concatenate(([0], dx, [0]))
        for i, xi in enumerate(x):
            yi = self.y[ np.logical_and( self.x>xi-dx[i], self.x<xi+dx[i+1] ) ]
            y[i] = yi.mean()
        self.x = x
        self.y = y

    def graficar(self, figax=None, **kwargs):
        """
        Graficar una curva.
        Args:
            figax: tuple with (fig, ax)
        Result:

        """
        if figax:
            fig,ax = figax
        else:
            fig, ax = plt.subplots(figsize=(8,6))
        kwargs['label'] = kwargs.get('label', self.nombre)
        ax.plot(self.x, self.y, **kwargs)
        return fig, ax

    def calcular_derivada(self, ispread=1):
        xd = .5 * (self.x[ispread:] + self.x[:-ispread])
        yd = ( self.y[ispread:] - self.y[:-ispread] ) / ( self.x[ispread:] - self.x[:-ispread] )
        self.d = self.__class__(self.nombre + '_d', xd, yd)

    def calcular_secante(self, v1, v2, axis=0):
        if axis==0: # v1 y v2 son valores en x
            x1 = v1
            x2 = v2
            y1 = np.interp(x1, self.x, self.y)
            y2 = np.interp(x2, self.x, self.y)
        elif axis==1:
            y1 = v1
            y2 = v2
            x1 = np.interp(y1, self.y, self.x)
            x2 = np.interp(y2, self.y, self.x)
        return (y2-y1) / (x2-x1)

class Curva_TvsD(Curva):
    def __init__(self, nombre, x, y):
        Curva.__init__(self, nombre, x, y)
        self.yield_defo = None
        self.yield_modtan = None
        self.yield_ten = None

    def descartar_rotura(self):
        max_ten = self.y.max()
        max_ten_index, = np.where( self.y==max_ten )[0]
        mask = self.x < self.x[max_ten_index]
        self.x = self.x[mask]
        self.y = self.y[mask]

    def calcular_yield(self, fr):
        # trabajo con la derivada de la curva tension-deformacion
        # si no tengo una derivada, la calculo
        if not self.d:
            self.calcular_derivada()
        # busco el maximo valor de modulo tangente
        y_max = self.d.y.max()
        # me fijo el indice donde ocurre
        y_max_index = np.where( self.d.y == y_max )[0][0]
        # selecciono la parte de la curva a la derecha de este valor
        mask = self.d.x >= self.d.x[y_max_index]
        x = self.d.x[mask]
        y = self.d.y[mask]
        # ahora encuentro el valor de x para el cual se cumple la condicion (y(x) = y_yield = fr*ymax)
        y_yield = y_max*fr
        x_yield = np.interp(y_yield, y[::-1], x[::-1])  # recorro hacia atras, sino no funciona interp
        # ahora me fijo los valores de deformacion y tension
        ten_yield = np.interp(x_yield, self.x, self.y)
        self.yield_defo = x_yield
        self.yield_modtan = y_yield
        self.yield_ten = ten_yield
        return x_yield, y_yield, ten_yield

    def get_modulo_por_defo(self, defo1, defo2):
        ten1, ten2 = np.interp([defo1, defo2], self.x, self.y)
        return (ten2-ten1) / (defo2-defo1)

    def get_modulo_por_ten(self, ten1, ten2):
        defo1, defo2 = np.interp([ten1, ten2], self.y, self.x)
        return (ten2-ten1) / (defo2-defo1)

    def calcular_modulos(self, ten1, ten2, defo1, defo2):
        self.Ee = self.get_modulo_por_ten(ten1, ten2)
        self.Epy = self.get_modulo_por_defo(defo1, defo2)

    def graficar(self, figax=None,
            yield_kwargs={'marker':'.', 'mec':'k', 'mfc':'w', 'ms':8},
            **kwargs):
        fig, ax = Curva.graficar(self, figax=figax, **kwargs)
        if self.yield_defo:
            ax.plot(self.yield_defo, self.yield_ten, lw=0, **yield_kwargs)
        return fig, ax

def round_up(x, d):
    aux = 10**d
    return np.ceil(x*aux) / float(aux)

class Conjunto(object):
    def __init__(self, nombre, c=[], mean=None, std=None, rstd=None):
        self.nombre = nombre
        self.c = c # lista de curvas
        self.mean = mean # curva promedio
        self.std = std # curva de desviacion estandar
        self.rstd = rstd # curva de desviacion estandar relativa
        self.reduced = False

    def reducir_por_x(self, x):
        for c in self.c:
            c.reducir_segun_x_promediando(x)
        self.reduced = True

    def reducir_por_dx(self, dx, d=None):
        x0max = max( [c.x[0] for c in self.c] )
        x1min = min( [c.x[-1] for c in self.c] )
        x0 = round_up(x0max, d) if d else x0max
        x = np.arange(x0, x1min, dx)
        for c in self.c:
            c.reducir_segun_x_promediando(x)
        self.reduced = True

    def armar_df(self):
        if not self.reduced:
            raise ValueError('Lista de curvas no reducida ni regularizada, imposible armar df!')
        df = pd.DataFrame()
        lens = [len(c.x) for c in self.c]
        x = self.c[lens.index(max(lens))].x
        df = pd.DataFrame(data={'x':x})
        for c in self.c:
            df_c = pd.DataFrame(data={c.nombre+'.y':c.y})
            df = pd.concat([df, df_c], axis=1)
        self.df = df

    def calcular_mean_std(self):
        cols = [col for col in self.df.columns if col[-2:]=='.y']
        nrows = self.df.count().min()
        self.df['mean'] = self.df.loc[:nrows-1,cols].mean(axis=1)
        self.df['std'] = self.df.loc[:nrows-1,cols].std(axis=1)
        self.df['rstd'] = self.df['std'] / self.df['mean']
        mean = self.df['mean'].values
        mask = ~np.isnan(mean)
        mean = mean[mask]
        std = self.df['std'].values[mask]
        x = self.df['x'].values[mask]
        rstd = self.df['rstd'].values[mask]
        self.mean = Curva_TvsD('mean', x, mean)
        self.std = Curva_TvsD('std', x, std)
        self.rstd = Curva_TvsD('rstd', x, rstd)

    def graficar(self, figax=None, **kwargs):
        if figax:
            fig, ax = figax
        else:
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
        for c in self.c:
            c.graficar(figax=(fig,ax), **kwargs)
            # ax.plot(c.x, c.y, label=c.nombre, **kwargs)
        return fig, ax

    def graficar_mean(self, figax=None, std=False, rstd=False, rstd_kwargs={}, **kwargs):
        if figax:
            fig, ax = figax
        else:
            fig, ax = plt.subplots(figsize=(8,6))
        self.mean.graficar(figax=(fig,ax), label=self.nombre, **kwargs)
        if std:
            ax.fill_between(self.mean.x, self.mean.y+self.std.y, self.mean.y-self.std.y,
                    alpha=0.2, zorder=1, label=None)
        if rstd:
            rstd.plot(self.rstd.x, self.rstd.y, **rstd_kwargs)
        return fig, ax

    def graficar_rstd(self, figax=None, **kwargs):
        if figax:
            fig, ax = figax
        else:
            fig, ax = plt.subplots(figsize=(8,6))
        #
        self.rstd.graficar(figax=(fig,ax), label=self.nombre, **kwargs)

class Conjunto_TvsD(Conjunto):
    def descartar_roturas(self):
        for c in self.c:
            c.descartar_rotura()

    def calcular_yields(self, fr):
        yield_ten_l = []
        for c in self.c:
            c.calcular_yield(fr)
            yield_ten_l.append(c.yield_ten)
        self.yield_ten = np.array(yield_ten_l)
        if self.mean:
            self.mean.calcular_yield(fr)

    def calcular_modulos(self, ten1, ten2, defo1, defo2):
        Ee_l = []
        Epy_l = []
        for c in self.c:
            c.calcular_modulos(ten1, ten2, defo1, defo2)
            Ee_l.append(c.Ee)
            Epy_l.append(c.Epy)
        self.Ee = np.array(Ee_l)
        self.Epy = np.array(Epy_l)
        if self.mean:
            self.mean.calcular_modulos(ten1, ten2, defo1, defo2)

class Conjuntos(object):
    def __init__(self, d={}):
        self.d = d

    def add_conjunto(self, c):
        self.d[c.nombre] = c

    def graficar_mean(self, figax=None, nombres=[], std=True, **kwargs):
        if figax:
            fig, ax = figax
        else:
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
        #
        nombres = nombres if nombres else self.d.keys()
        for nombre in nombres:
            c = self.d[nombre]
            c.graficar_mean(figax=(fig,ax), std=std, **kwargs)
        #
        return fig, ax

    def graficar_rstd(self, figax=None, nombres=[], **kwargs):
        if figax:
            fig, ax = figax
        else:
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
        #
        nombres = nombres if nombres else self.d.keys()
        for nombre in nombres:
            c = self.d[nombre]
            c.graficar_rstd(figax=(fig,ax))
        #
        return fig, ax