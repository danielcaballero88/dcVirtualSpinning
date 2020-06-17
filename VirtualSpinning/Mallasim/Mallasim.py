# Global imports
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
# Local imports
from .Fibras import Fibras
from .Nodos import Nodos
from VirtualSpinning.aux import find_string_in_file

# Parameters
DELTA = 1.e-4
DELTA21 = 1. / (2. * DELTA)
DELTAX = DELTA * np.array([1., 0.], dtype=float)
DELTAY = DELTA * np.array([0., 1.], dtype=float)


class Mallasim(object):
    """
    Mallasim es una clase simplificada de mallas para tratar
    las mallas desde un punto de vista mecanico, no geometrico
    Aqui cada fibra no es el cuerpo depositado sino
    el cuerpo comprendido entre dos puntos de union (interseccion)
    """

    def __init__(self, L, Dm, nodos, fibras):
        self.L = L
        self.Dm = Dm
        self.nodos = nodos
        self.fibras = fibras
        self.status_deformed = False
        self.Fmacro = None
        self.Tmacro = None
        self.Ncapas = None

    @classmethod
    def leer_desde_archivo(cls, nomarchivo):
        fid = open(nomarchivo, 'r')
        # primero leo los parametros
        target = "*parametros"
        ierr = find_string_in_file(fid, target, True)
        L = float(next(fid))
        Dm = float(next(fid))
        Nc = int(next(fid))
        nparam = int(next(fid))
        svals = next(fid).split()
        param_in = [float(val) for val in svals]
        # me fijo su es una malla deformada
        target = "*deformacion"
        ierr = find_string_in_file(fid, target, mandatory=False)
        status_deformed = False
        if ierr == 0:
            status_deformed = True
            linea = next(fid)
            F11, F21, F12, F22 = [float(item) for item in linea.split()]
            Fmacro = np.array([[F11, F12], [F21, F22]])
            linea = next(fid)
            T11, T21, T12, T22 = [float(item) for item in linea.split()]
            Tmacro = np.array([[T11, T12], [T21, T22]])
        else:
            status_deformed = False
            F11, F21, F12, F22 = 1.0, 0.0, 0.0, 1.0
            Fmacro = np.array([[F11, F12], [F21, F22]])
            Tmacro = np.array([[0.0, 0.0], [0.0, 0.0]])
        # luego busco coordenadas
        target = "*coordenadas"
        ierr = find_string_in_file(fid, target, True)
        num_r = int(next(fid))
        coors0 = list()
        coors = list()
        tipos = list()
        # -
        L_2 = 0.5 * L
        if status_deformed:
            L_2_def_x = 0.5 * L * F11
            L_2_def_y = 0.5 * L * F22
        else:
            L_2_def_x = L_2
            L_2_def_y = L_2
        # -
        for i in range(num_r):
            j, t, x0, y0, x, y = (float(val) for val in next(fid).split())
            x0 = x0 - L_2
            y0 = y0 - L_2
            if status_deformed:
                x = x - L_2_def_x
                y = y - L_2_def_y
            else:
                x = x0
                y = y0
            tipos.append(int(t))
            coors0.append([x0, y0])
            coors.append([x, y])
        # luego las fibras
        target = "*fibras"
        ierr = find_string_in_file(fid, target, True)
        num_f = int(next(fid))
        fibs = list()
        ds = list()
        letes0 = list()
        lamsr = list()
        lamps = list()
        brokens = list()
        for i in range(num_f):
            svals = next(fid).split()
            j = int(svals[0])
            d = float(svals[1])
            lete0 = float(svals[2])
            lamr0 = float(svals[3])
            lamp = float(svals[4])
            broken = False if svals[5] == "F" else True
            n0 = int(svals[6])
            n1 = int(svals[7])
            fibs.append([n0, n1])
            ds.append(d)
            letes0.append(lete0)
            lamsr.append(lamr0)
            lamps.append(lamp)
            brokens.append(broken)
        fid.close()
        # ---
        # ahora coloco las variables en mi objeto malla simplificada
        # nodos con coordenadas y tipos
        nodos = Nodos(num_r, coors0, tipos)
        nodos.r = np.array(coors, dtype=float)
        # subfibras
        letes0 = np.array(letes0, dtype=float)
        lamsr = np.array(lamsr, dtype=float)
        locos = letes0 * lamsr
        param = np.zeros((num_f, nparam), dtype=float)
        param[:, 0:] = param_in
        fibras = Fibras(num_f, fibs, ds, letes0, lamsr, lamps, brokens, param)
        # mallita
        m = cls(L, Dm, nodos, fibras)
        m.Ncapas = Nc
        if status_deformed:
            m.status_deformed = status_deformed
            m.Fmacro = Fmacro
            m.Tmacro = Tmacro
        return m

    def escribir_en_archivo(self, nomarchivo):
        fid = open(nomarchivo, "w")
        # primero escribo L, dl y dtheta
        dString = "*Parametros \n"
        fmt = "{:20.8e}"
        dString += fmt.format(self.L) + "\n"
        dString += fmt.format(self.Dm) + "\n"
        dString += "{:3d}".format(self.Ncapas) + "\n"
        nparam = np.shape(self.fibras.param)[1]
        dString += "{:3d}".format(nparam) + "\n"
        dString += "".join(fmt.format(val) for val in self.fibras.param[0]) + "\n"
        fid.write(dString)
        # escribo los nodos: indice, tipo, y coordenadas
        dString = "*Coordenadas \n" + str(self.nodos.n) + "\n"
        fid.write(dString)
        for n in range(self.nodos.n):
            dString = "{:12d}".format(n)
            dString += "{:2d}".format(self.nodos.t[n])
            dString += "".join("{:20.8e}".format(val) for val in self.nodos.r0[n]) + "\n"
            fid.write(dString)
        # ---
        # sigo con las fibras: indice, dl, d, dtheta, y segmentos (conectividad)
        dString = "*Fibras \n" + str(self.fibras.n) + "\n"
        fid.write(dString)
        for f, (n0, n1) in enumerate(self.fibras.con):
            dString = "{:12d}".format(f)  # indice
            dString += "{:20.8e}{:20.8e}{:20.8e}".format(self.fibras.ds[f], self.fibras.letes0[f],
                                                         self.fibras.lamsr[f])  # d, lete y lamr
            dString += "{:12d}{:12d}".format(n0, n1) + "\n"
            fid.write(dString)
        # ---
        fid.close()

    def calcular_fuerzas(self, longout=False):
        return self.fibras.calcular_fuerzas(self.nodos.r, longout)

    def deformar_frontera(self, F):
        self.nodos.r[self.nodos.mf] = np.matmul(self.nodos.r0[self.nodos.mf], np.transpose(F))

    def calcular_dr(self):
        dr = np.zeros(np.shape(self.nodos.r), dtype=float)
        r1 = self.nodos.r.copy()  # posiciones actualizadas en iteraciones
        for i in range(10):
            A, b = self.calcular_A_b(r1)
            dr = np.linalg.solve(A, b)
            residuos = np.matmul(A, dr) - b
            error = np.max(np.abs(np.sqrt(np.sum(dr * dr, axis=1))))
            print(i, ":", error)
            r1 += dr.reshape(-1, 2)
            if error < 1.e-10:
                break
        return r1

    def calcular_A_b(self, r1):
        """ calcula la matriz tangente A y el vector de cargas b
        alrededor de las posiciones dadas por r1
        que es tipicamente el array de posiciones en la iteracion previa """
        # seteo matrices
        nG = self.nodos.n * 2
        matG = np.zeros((nG, nG), dtype=float)
        vecG = np.zeros((nG, 1), dtype=float)
        nL = 2  # tengo que armar la matriz tangente respecto de un solo nodo (hay doble simetria)
        matL = np.zeros((nL, nL), dtype=float)
        vecL = np.zeros((nL, 1), dtype=float)
        for f, (n0, n1) in enumerate(self.fibras.con):
            if f == 2:
                pass
            r_n0 = r1[n0]
            r_n1 = r1[n1]
            r_n0_px = r1[n0] + DELTAX
            r_n0_mx = r1[n0] - DELTAX
            r_n0_py = r1[n0] + DELTAY
            r_n0_my = r1[n0] - DELTAY
            F_c = self.fibras.calcular_fuerza(f, r_n0, r_n1)
            F_mx = self.fibras.calcular_fuerza(f, r_n0_mx, r_n1)
            F_px = self.fibras.calcular_fuerza(f, r_n0_px, r_n1)
            F_my = self.fibras.calcular_fuerza(f, r_n0_my, r_n1)
            F_py = self.fibras.calcular_fuerza(f, r_n0_py, r_n1)
            dFdx = (F_px - F_mx) * DELTA21
            dFdy = (F_py - F_my) * DELTA21
            matL[:, 0] = dFdx[:, 0]
            matL[:, 1] = dFdy[:, 0]
            vecL = - F_c
            # ahora a ensamblar
            # primero el vector de cargas
            row = n0 * 2
            col = n1 * 2
            vecG[row:row + 2] += vecL
            vecG[col:col + 2] += -vecL
            # luego matriz local va a 4 submatrices de la global
            # primero en el nodo 0
            row = n0 * 2
            col = n0 * 2
            matG[row:row + 2, col:col + 2] += matL
            # luego lo mismo en el nodo 1
            row = n1 * 2
            col = n1 * 2
            matG[row:row + 2, col:col + 2] += matL
            # luego las cruzadas
            row = n0 * 2
            col = n1 * 2
            matG[row:row + 2, col:col + 2] += - matL
            row = n1 * 2
            col = n0 * 2
            matG[row:row + 2, col:col + 2] += - matL
        # ahora las condiciones de dirichlet
        for n in range(self.nodos.n):
            if self.nodos.t[n] == 1:
                ix = 2 * n
                iy = 2 * n + 1
                matG[ix, :] = 0.
                matG[ix, ix] = 1.
                vecG[ix] = 0.
                matG[iy, :] = 0.
                matG[iy, iy] = 1.
                vecG[iy] = 0.
        # fin
        return matG, vecG

    def pre_graficar_bordes(self, fig, ax, byn=False):
        # deformacion afin del borde
        r_corners = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5], [-.5, -.5]]) * self.L
        if self.status_deformed:
            Fmacro = self.Fmacro
            r_corners = np.matmul(r_corners, np.transpose(Fmacro))
        # seteo limites
        lim_left = np.min(r_corners[:, 0])
        lim_right = np.max(r_corners[:, 0])
        lim_bottom = np.min(r_corners[:, 1])
        lim_top = np.max(r_corners[:, 1])
        margen = 0.1 * self.L
        ax.set_xlim(left=lim_left - margen, right=lim_right + margen)
        ax.set_ylim(bottom=lim_bottom - margen, top=lim_top + margen)
        # dibujo los bordes del rve
        plt_fron = ax.plot(r_corners[:, 0], r_corners[:, 1], linestyle=":", c="gray")

    def pre_graficar_0(self, fig, ax, lamr_min=None, lamr_max=None, plotnodos=False, maxnfibs=500, colorbar=False):
        mi_cm = plt.cm.jet
        lamsr = self.fibras.lamsr
        if lamr_min is None:
            lamr_min = np.min(lamsr)
        if lamr_max is None:
            lamr_max = np.max(lamsr)
        sm = plt.cm.ScalarMappable(cmap=mi_cm, norm=plt.Normalize(vmin=lamr_min, vmax=lamr_max))
        nfibs = self.fibras.n
        if maxnfibs < nfibs:
            nfibs = maxnfibs
        for f in range(nfibs):
            n0, n1 = self.fibras.con[f]
            # linea inicial
            x0, y0 = self.nodos.r0[n0]
            x1, y1 = self.nodos.r0[n1]
            c = sm.to_rgba(lamsr[f])
            ax.plot([x0, x1], [y0, y1], ls="-", c=c)
        if plotnodos:
            xnods = self.nodos.r0[:, 0]
            ynods = self.nodos.r0[:, 1]
            ax.plot(xnods, ynods, linewidth=0, marker="o", mec="k", mfc="w", markersize=4)
        if colorbar:
            sm._A = []
            fig.colorbar(sm)

    @staticmethod
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    @staticmethod
    def calcular_tensiones_fibras(lams, Et, Eb, lamrs, lamps):
        nfibs = len(lams)
        tens = np.zeros(nfibs, dtype=float)
        for f in range(nfibs):
            lamrL = lamrs[f] * lamps[f]
            if lams[f] <= lamrL:
                tens[f] = Eb * (lams[f] / lamps[f] - 1.)
            else:
                tens[f] = Eb * (lamrL - 1.) + Et * (lams[f] / lamrL - 1.)
        return tens

    def pre_graficar(self, fig, ax, linewidth=2,
                     lam_min=None, lam_max=None,
                     maxnfibs=5000, byn=False, color_por="nada", barracolor=False, colormap="jet", colores_cm=None,
                     ncolores_cm=100,
                     afin=True, colorafin="gray", linewidthafin=2):
        print("pregraficando mallita")
        drs, longs, lams = self.fibras.calcular_drs_letes_lams(self.nodos.r)
        lamsr = self.fibras.lamsr
        lamps = self.fibras.lamps
        lams_ef = lams[:, 0] / lamsr / lamps
        if self.status_deformed:
            Fmacro = self.Fmacro
        else:
            Fmacro = np.array([[1., 0.], [0., 1.]])
        if byn:
            mi_cm = plt.cm.gray_r
            # lo trunco para que no incluya el cero (blanco puro que no hace contraste con el fondo)
            mi_cm = self.truncate_colormap(mi_cm, 0.4, 1.0)
        else:
            if colores_cm is not None:
                mi_cm = colors.LinearSegmentedColormap.from_list("mi_cm", colores_cm, N=ncolores_cm)
            elif colormap == "jet":
                mi_cm = plt.cm.jet
            elif colormap == "rainbow":
                mi_cm = plt.cm.rainbow
            elif colormap == "prism":
                mi_cm = plt.cm.prism
            elif colormap == "Dark2":
                mi_cm = plt.cm.Dark2
        if color_por == "lamr":
            if lam_min is None:
                lam_min = np.min(lamsr)
            if lam_max is None:
                lam_max = np.max(lamsr)
            sm = plt.cm.ScalarMappable(cmap=mi_cm, norm=plt.Normalize(vmin=lam_min, vmax=lam_max))
        elif color_por == "lam":
            if lam_min is None:
                lam_min = np.min(lams)
            if lam_max is None:
                lam_max = np.max(lams)
            sm = plt.cm.ScalarMappable(cmap=mi_cm, norm=plt.Normalize(vmin=lam_min, vmax=lam_max))
        elif color_por == "lam_ef":
            if lam_min is None:
                lam_min = np.min(lams_ef)
            if lam_max is None:
                lam_max = np.max(lams_ef)
            sm = plt.cm.ScalarMappable(cmap=mi_cm, norm=plt.Normalize(vmin=lam_min, vmax=lam_max))
        elif color_por == "tension":
            tensiones = self.calcular_tensiones_fibras(lams, 2.9e3, 2.0, lamsr, lamps)
            if lam_min is None:
                lam_min = np.min(tensiones)
            if lam_max is None:
                lam_max = np.max(tensiones)
            sm = plt.cm.ScalarMappable(cmap=mi_cm, norm=plt.Normalize(vmin=lam_min, vmax=lam_max))
        elif color_por == "fibra":
            sm = plt.cm.ScalarMappable(cmap=mi_cm, norm=plt.Normalize(vmin=0, vmax=self.fibras.n - 1))
        elif color_por == "reclutamiento":
            lams_ef = lams[:, 0] / lamsr / lamps
            sm = plt.cm.ScalarMappable(cmap=mi_cm, norm=plt.Normalize(vmin=0, vmax=self.fibras.n - 1))  # al pedo
        elif color_por == "nada":
            sm = plt.cm.ScalarMappable(cmap=mi_cm, norm=plt.Normalize(vmin=0, vmax=self.fibras.n - 1))  # al pedo

        if afin is False or Fmacro is None:
            # si se da Fmacro r0 se modifica para que sean las coordenadas de deformacion afin
            # de lo contrario quedan las iniciales
            r0 = self.nodos.r0
        else:
            r0 = np.matmul(self.nodos.r0, np.transpose(Fmacro))
        nfibs = self.fibras.n
        if maxnfibs < nfibs:
            nfibs = maxnfibs
        pctp = np.arange(0., 101., 10.).tolist()
        pcpd = np.zeros(len(pctp), dtype=bool).tolist()
        for f in range(nfibs):
            # imprimo el porcentaje de completado
            pc = round(float(f) / nfibs * 100., 0)
            if pc in pctp:
                ipc = pctp.index(pc)
                if not pcpd[ipc]:
                    print("{:4.0f}% ".format(pc), end='')
                    pcpd[ipc] = True
            # ---
            n0, n1 = self.fibras.con[f]
            if afin and Fmacro is not None:
                # linea inicial (afin si se da Fmacro)
                x0, y0 = r0[n0]
                x1, y1 = r0[n1]
                ax.plot([x0, x1], [y0, y1], ls=":", c=colorafin, linewidth=linewidthafin)
            # linea final
            x0, y0 = self.nodos.r[n0]
            x1, y1 = self.nodos.r[n1]
            if color_por == "lamr":
                c = sm.to_rgba(lamsr[f])
            elif color_por == "lam":
                c = sm.to_rgba(lams[f, 0])
            elif color_por == "lam_ef":
                c = sm.to_rgba(lams_ef[f])
            elif color_por == "tension":
                c = sm.to_rgba(tensiones[f])
            elif color_por == "fibra":
                c = sm.to_rgba(f)
            elif color_por == "reclutamiento":
                if lams_ef[f] > 1. + 1.e-6:
                    c = "red"
                else:
                    c = "blue"
            elif color_por == "nada":
                c = "k"
            # --
            # if self.fibras.brokens[f]:
            #     c = "gray"
            # elif lams_ef[f] < 1.00001:
            #     c = "k"
            # --
            ax.plot([x0, x1], [y0, y1], ls="-", linewidth=linewidth, c=c)
        print()
        if barracolor:
            sm._A = []
            fig.colorbar(sm)
