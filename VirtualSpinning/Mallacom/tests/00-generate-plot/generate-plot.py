from numpy import pi
from matplotlib import pyplot as plt
from pathlib import Path
import time
from VirtualSpinning.Mallacom.Mallacom import Mallacom

CWD = Path.cwd()
FILE = Path(__file__)
DIR = FILE.parent


def main():
    start = time.time()
    print('testing: Mallacom -> generate-plot')
    print('CWD: ', CWD)
    print('DIR: ', DIR)
    print('__file__: ', FILE)
    Dm = 1.
    L = 100. * Dm
    fundisor = None
    dl = 5. * Dm
    devang = 10. * pi / 180.
    volfrac = 0.1
    ncaps = 2
    mc = Mallacom(L, Dm, volfrac, dl, devang)
    for _i in range(1, ncaps + 1):
        mc.make_capa(dl, Dm, devang, volfrac, orient_distr=fundisor)
    archivo = DIR / 'malla_test_basic.txt'
    mc.guardar_en_archivo(archivo)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    mc.marco.graficar(fig, ax)
    mc.pre_graficar_fibras(fig, ax, color_por="capa", byn=True, linewidth=1.5)
    print('time: ', time.time() - start)
    plt.show()


if __name__ == '__main__':
    main()
