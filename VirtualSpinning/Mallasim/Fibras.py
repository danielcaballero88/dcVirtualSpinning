import numpy as np


class Fibras(object):
    def __init__(self, n, con, ds, letes0, lamsr, lamps, brokens, param):
        self.n = n
        self.con = np.array(con, dtype=int)
        self.ds = np.array(ds, dtype=float)  # diametros
        self.letes0 = np.array(letes0, dtype=float)
        self.lamsr = np.array(lamsr, dtype=float)
        self.lamps = np.array(lamps, dtype=float)
        self.brokens = np.array(brokens, dtype=bool)
        self.param = np.array(param, dtype=float)
        self.drs0 = None
        self.__mr = np.zeros((n, 1), dtype=bool)
        self.__me = np.zeros((n, 1), dtype=bool)
        self.__fzas = np.zeros((n, 1), dtype=float)
        self.__fzasv = np.zeros((n, 2), dtype=float)

    def calcular_drs0(self, r0):
        self.drs0 = r0[self.con[:, 1]] - r0[self.con[:, 0]]
        return self.drs0

    def calcular_fzasr(self):
        self.fzas_r = np.zeros((self.n, 2), dtype=float)
        self.fzas_r = self.param[:, 2, None] * (self.lamsr - 1.)

    def calcular_drs_letes_lams(self, r):
        drs = r[self.con[:, 1]] - r[self.con[:, 0]]
        longs = np.sqrt(np.sum(drs * drs, axis=1, keepdims=True))
        lams = longs / self.letes0[:, None]
        return drs, longs, lams

    def calcular_lams(self, r):
        *_, lams = self.calcular_drs_letes_lams(r) 
        return lams

    def calcular_fuerzas(self, r, longout=False):
        drs, longs, lams = self.calcular_drs_letes_lams(r)
        lams_r = self.lamsr
        ks1 = self.param[:, 1, None]
        ks2 = self.param[:, 2, None]
        self.__mr = np.greater(lams, lams_r)  # mask rectas
        self.__me = np.logical_not(self.__mr)  # mask enruladas
        i = self.__mr
        self.__fzas[i] = self.fzas_r[i] + ks1[i] * (lams[i] / lams_r[i] - 1.)
        i = self.__me
        self.__fzas[i] = ks2[i] * (lams[i] - 1.)
        self.__fzasv = self.__fzas / longs * drs
        if not longout:
            return self.__fzasv
        else:
            return drs, longs, lams, self.__fzas, self.__fzasv

    def calcular_fuerza(self, f, r_n0, r_n1, longout=False):
        dr = r_n1 - r_n0
        long = np.sqrt(np.sum(dr * dr))
        lam = long / self.letes0[f]
        lam_r = self.lamsr[f]
        k1, k2 = self.param[f]
        fza_r = k2 * (lam_r - 1.)
        if lam <= lam_r:
            fza = k2 * (lam - 1.)
        else:
            fza = fza_r + k1 * (lam / lam_r - 1.)
        fzav = fza / long * dr[:, None]
        if not longout:
            return fzav
        else:
            return dr, long, lam, fza, fzav