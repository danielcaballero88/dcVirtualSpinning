# Global imports
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import cm
# Local imports
from .Fibras import Fibras
from .Nodos import Nodos
from .Marco import Marco
from VirtualSpinning.aux import find_string_in_file

# Parameters
DELTA = 1.e-4
DELTA21 = 1. / (2. * DELTA)
DELTAX = DELTA * np.array([1., 0.], dtype=float)
DELTAY = DELTA * np.array([0., 1.], dtype=float)


SMALL_SIZE = 12
MEDIUM_SIZE = 20
BIGGER_SIZE = 24
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


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
        self.marco = Marco(L)
        self.nparcon = 0 
        self.parcon = []

    @classmethod
    def leer_de_archivo(cls, nomarchivo):
        fid = open(nomarchivo, 'r')
        # primero leo los parametros
        target = "*parametros"
        ierr = find_string_in_file(fid, target, True)
        L = float(next(fid))
        Dm = float(next(fid))
        Nc = int(next(fid))
        # me fijo su es una malla deformada
        target = "*deformacion"
        ierr = find_string_in_file(fid, target, mandatory=False)
        status_deformed = False
        if ierr == 0:
            # si la malla esta deformada leo su Fmacro y Tmacro
            status_deformed = True
            linea = next(fid)
            F11, F21, F12, F22 = [float(item) for item in linea.split()]
            Fmacro = np.array([[F11, F12], [F21, F22]])
            linea = next(fid)
            T11, T21, T12, T22 = [float(item) for item in linea.split()]
            Tmacro = np.array([[T11, T12], [T21, T22]])
            # tambien sus parametros constitutivos
            nparcon = int(next(fid))
            svals = next(fid).split()
            parcon = [float(val) for val in svals]
        else:
            # si no esta deformada pongo valores default
            status_deformed = False
            F11, F21, F12, F22 = 1.0, 0.0, 0.0, 1.0
            Fmacro = np.array([[F11, F12], [F21, F22]])
            Tmacro = np.array([[0.0, 0.0], [0.0, 0.0]])
            nparcon = 0 
            parcon = []
        # luego busco coordenadas
        target = "*coordenadas"
        ierr = find_string_in_file(fid, target, True)
        num_r = int(next(fid))
        coors0 = list()
        coors = list()
        tipos = list()
        for _i in range(num_r):
            _j, t, x0, y0, x, y = (float(val) for val in next(fid).split())
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
        for _i in range(num_f):
            svals = next(fid).split()
            _j = int(svals[0])
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
        _locos = letes0 * lamsr
        param = np.zeros((num_f, nparcon), dtype=float)
        param[:, 0:] = parcon
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

    def pre_graficar_bordes(self, fig, ax, limites={}):
        # deformacion afin del borde
        if self.status_deformed:
            self.marco.deformar(self.Fmacro)
        # graficar
        self.marco.graficar(fig, ax, limites)

    def pre_graficar_0(self, fig, ax, lamr_min=None, lamr_max=None, plotnodos=False, maxnfibs=500, colorbar=False):
        mi_cm = cm.get_cmap('jet')
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
    def trunc_cm(cmap, minval=0.0, maxval=1.0, n=100):
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

    def pre_graficar(self, fig, ax,
                    maxnfibs=5000,
                    cby='nada', cvmin=None, cvmax=None,
                    byn=False, cmap='jet', cbar=True, cbar_ticks=None,
                    clist=None, N_clist=20,
                    plot_inic=False, parplotinic={},
                    plot_afin=False, parplotafin={},
                    plot_broken=False, plot_enrul=False,
                    parplot={}):
        """
        Metodo para graficar el rve de una Mallasim
        """

        print("Pregraficando Mallasim")

        # Preparo el mapa de colores
        if byn:
            mi_cm = cm.get_cmap('gray_r') # 0=blanco, 1=negro
            # Lo trunco para que no incluya el cero 
            # (blanco puro que no hace contraste con el fondo)
            mi_cm = self.trunc_cm(mi_cm, 0.4, 1.0) 
        elif clist is not None:
            mi_cm = colors.LinearSegmentedColormap.from_list("mi_cm", clist, N=N_clist)
        else:
            mi_cm = cm.get_cmap(cmap)  # default: jet

        # Calculo variables de las fibras
        lams = self.fibras.calcular_lams(self.nodos.r)[:,0]
        lamsr = self.fibras.lamsr
        lamps = self.fibras.lamps
        lams_ef = lams / lamsr / lamps

        # Me fijo segun que variable coloreo
        if cby == 'lamr':
            cvar = lamsr
        elif cby == 'lam':
            cvar = lams
        elif cby == 'lam_ef':
            cvar = lams_ef
        elif cby == 'lamp':
            cvar = lamps
        elif cby == 'fibra':
            cvar = list(range(self.fibras.n))
            if cvmin is None: cvmin = 0
            if cvmax is None: cvmax = self.fibras.n - 1
        elif cby == 'reclutamiento': 
            mask = lams_ef > 1.0 + 1.0e-6
            cvar = np.zeros(self.fibras.n, dtype=float)
            cvar[mask] = 1.
            cvmin = 0
            cvmax = 1
        elif cby == 'nada':
            cvar = np.zeros(self.fibras.n, dtype=float)
            cvmin = 0
            cvmax = 1
            mi_cm = cm.get_cmap('gray')

        # Pongo los min y max defaults en cayo de que no esten seteados ya
        if cvmin is None: cvmin = cvar.min()
        if cvmax is None: cvmax = cvar.max()

        # Armo un ScalarMappable Colormap para poder colorear segun variable
        sm = plt.cm.ScalarMappable(cmap=mi_cm, norm=plt.Normalize(vmin=cvmin, vmax=cvmax))

        # Decido si ademas de la malla deformada grafico la deformacion afin para comparar
        if plot_afin: rafin = np.matmul(self.nodos.r0, np.transpose(self.Fmacro))
        # Idem para la configuracion inicial 
        if plot_inic: rinic = self.nodos.r0

        # Me fijo si tengo un maximo numero de fibras a graficar
        nfibs = self.fibras.n
        if maxnfibs < nfibs: nfibs = maxnfibs

        # Armo algunas variables para medir porcentaje de avance
        pctp = np.arange(0., 101., 10.).tolist()
        pcpd = np.zeros(len(pctp), dtype=bool).tolist()

        # Empiezo el bucle para graficar
        for f in range(nfibs):

            # Imprimo el porcentaje de completado
            pc = round(float(f) / nfibs * 100., 0)
            if pc in pctp:
                ipc = pctp.index(pc)
                if not pcpd[ipc]:
                    print("{:4.0f}% ".format(pc), end='')
                    pcpd[ipc] = True
            
            # Grafico la fibra f
            n0, n1 = self.fibras.con[f]

            # Grafico configuracion inicial (para comparar)
            if plot_inic: 
                x0, y0 = rinic[n0]
                x1, y1 = rinic[n1]
                ax.plot([x0, x1], [y0, y1], **parplotinic)

            # Grafico configuracion afin (para comparar)
            if plot_afin:
                x0, y0 = rafin[n0]
                x1, y1 = rafin[n1]
                # algunos defaults
                parplotafin['c'] = parplotafin.get('c', 'gray')
                parplotafin['ls'] = parplotafin.get('ls', ':')
                ax.plot([x0, x1], [y0, y1], **parplotafin)
            
            # Grafico configuracion deformada
            x0, y0 = self.nodos.r[n0]
            x1, y1 = self.nodos.r[n1]
            c = sm.to_rgba(cvar[f])
            # --
            # Si la fibra se ha roto o aun esta sin reclutar le pongo otro color
            if plot_broken:
                if self.fibras.brokens[f]:
                    c = "gray"
            if plot_enrul:
                if lams_ef[f] < 1.00001:
                    c = "k"
            # --
            # La agrego al plot
            ax.plot([x0, x1], [y0, y1], c=c, **parplot)
        print()

        # Agrego la colorbar si asi lo quise
        if cbar:
            sm._A = []
            fig.colorbar(sm)

        # Fin
