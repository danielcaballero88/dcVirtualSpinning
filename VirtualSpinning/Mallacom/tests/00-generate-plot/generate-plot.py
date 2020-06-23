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
    archivo = DIR / 'temp' / f'{mc.name}.txt'
    mc.guardar_en_archivo(archivo)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    mc.marco.graficar(fig, ax)
    mc.pre_graficar_fibras(fig, ax, color_por="capa", byn=True, linewidth=1.5)
    print('time: ', time.time() - start)
    fig.savefig(DIR / 'temp' / 'malla.png')

    fig, ax = plt.subplots(figsize=(8,6))
    mc = Mallacom.leer_de_archivo(archivo) 
    mc.marco.graficar(fig, ax) 
    mc.pre_graficar_fibras(fig, ax, color_por="capa", byn=True, linewidth=1.5)

    plt.show()


if __name__ == '__main__':
    main()
