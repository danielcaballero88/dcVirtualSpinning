import numpy as np
import scipy.stats as stats


class Random(object):
    def __init__(self, lower=0., upper=1.):
        self.lower = lower
        self.upper = upper
        self.range = upper - lower

    def __call__(self):
        return self.lower + np.random.rand() * self.range


class NormalTruncada(object):
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


class ValoresPreestablecidos(object):
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
