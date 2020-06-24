from numpy import pi
from matplotlib import pyplot as plt
from pathlib import Path
import time
from VirtualSpinning.Mallacom.Mallacom import Mallacom
from VirtualSpinning import odf


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
    fdo = odf.NormTr(n=100, loc=.5*pi, scale=.1*pi, lower=0., upper=pi)
    ls = 5. * D
    dth = 10. * pi / 180.
    vf = 0.1
    nc = 2
    params = {
        'L': L,
        'D': D,
        'vf': vf,
        'ls': ls,
        'dth': dth,
        'nc': nc,
        'fdo': fdo,
        'nm': 1
    }
    mc = Mallacom(**params, name='malla')
    mc.make_malla()
    archivo = DIR / f'{mc.name}.txt'
    mc.guardar_en_archivo(archivo)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    mc.marco.graficar(fig, ax)
    mc.pre_graficar_fibras(fig, ax, cby="capa", byn=True)
    print('time: ', time.time() - start)
    fig.savefig(DIR / 'malla.png')
    plt.show()


if __name__ == '__main__':
    main()
