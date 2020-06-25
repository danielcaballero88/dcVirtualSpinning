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


def main():
    start = time.time()
    print('testing: Mallacom -> generate-plot')
    print('CWD: ', CWD)
    print('DIR: ', DIR)
    print('__file__: ', FILE)

    # Armo los parametros como listas a combinar
    params = {
        'D' : [1.],
        'L' : np.array([50, 100.]),
        'fdo' : [None],
        'ls' : np.array([5.]),
        'dth' : np.array([10.]) * PI / 180.,
        'vf' : [0.1],
        'nc' : [2],
        'nm' : list(range(1, 2+1))    
    }
    # Convierto a valores dimensionales
    params['L'] *= params['D']
    params['ls'] *= params['D']

    # Armo las mallas
    mcs = Mallacom.make_combinaciones(params)
    print(time.time() - start)

    # Las grafico
    for mc in mcs:
        print(mc.param)
        archivo = DIR / f'{mc.name}.txt'
        mc.guardar_en_archivo(archivo)
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        mc.marco.graficar(fig, ax)
        mc.pre_graficar_fibras(fig, ax, cby="capa", byn=True)
        fig.savefig(DIR / f'{mc.name}.png')

    print('time: ', time.time() - start)
    plt.show()


if __name__ == '__main__':
    main()
