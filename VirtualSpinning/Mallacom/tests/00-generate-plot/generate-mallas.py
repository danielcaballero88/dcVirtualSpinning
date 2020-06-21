import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import time
from VirtualSpinning.Mallacom.Mallacom import Mallacom
from VirtualSpinning import odf


CWD = Path.cwd()
FILE = Path(__file__)
DIR = FILE.parent
PI = np.pi

def assemble_params(l_L, l_dl, l_dth, l_nc, nm):
    """
    A partir de listas de parametros (si o si listas) armo una
    lista con los diccionarios para armar las mallas, 
    y tambien los nombres de las mallas
    """
    for L in l_L:
        for dl in l_dl:
            for dth in l_dth:
                for nc in l_nc:
                    for im in range(1,nm+1):
                        print(f'nc={nc:05d}  L={L:08.2f}  dth={dth:05.2f}  dl={dl:05.2f}  im={im:03d}')



def main():
    start = time.time()
    print('testing: Mallacom -> generate-plot')
    print('CWD: ', CWD)
    print('DIR: ', DIR)
    print('__file__: ', FILE)
    Dm = 1.
    L = np.array([50, 100.]) * Dm
    fdo = odf.NormTr(n=100, loc=.5*PI, scale=.1*PI, lower=0., upper=PI)
    dl = np.array([5.]) * Dm
    devang = np.array([10.]) * PI / 180.
    volfrac = 0.1
    ncaps = [2]
    nm = 2
    assemble_params(L, dl, devang, ncaps, nm)
    mc = Mallacom(L, Dm, volfrac, dl, devang)
    for _i in range(1, ncaps + 1):
        mc.make_capa(dl, Dm, devang, volfrac, orient_distr=fdo)
    archivo = DIR / 'malla_generate-aligned.txt'
    mc.guardar_en_archivo(archivo)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    mc.marco.graficar(fig, ax)
    mc.pre_graficar_fibras(fig, ax, color_por="capa", byn=True, linewidth=1.5)
    print('time: ', time.time() - start)
    plt.show()


if __name__ == '__main__':
    main()
