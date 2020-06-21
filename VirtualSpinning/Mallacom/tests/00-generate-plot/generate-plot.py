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
    D = 1.
    L = 100. * D
    fdo = None
    ls = 5. * D
    dth = 10. * pi / 180.
    vf = 0.1
    nc = 2
    mc = Mallacom(L, D, vf, ls, dth, nc)
    for _i in range(1, nc + 1):
        mc.make_capa(fdo=fdo)
    archivo = DIR / 'temp' / 'malla_generate-plot.txt'
    mc.guardar_en_archivo(archivo)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    mc.marco.graficar(fig, ax)
    mc.pre_graficar_fibras(fig, ax, color_por="capa", byn=True, linewidth=1.5)
    print('time: ', time.time() - start)
    plt.show()


if __name__ == '__main__':
    main()
