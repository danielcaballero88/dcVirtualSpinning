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


def assemble_params(l_L, l_ls, l_dth, l_nc, l_nm):
    """
    A partir de listas de parametros (si o si listas) armo una
    lista con los diccionarios para armar las mallas, 
    y tambien los nombres de las mallas
    """
    params = []
    for L in l_L:
        for ls in l_ls:
            for dth in l_dth:
                for nc in l_nc:
                    for nm in l_nm:
                        d = {
                            'L':L,
                            'ls': ls,
                            'dth': dth,
                            'nc': nc,
                            'nm': nm
                            }
                        params.append(d)
    return params


def main():
    start = time.time()
    print('testing: Mallacom -> generate-plot')
    print('CWD: ', CWD)
    print('DIR: ', DIR)
    print('__file__: ', FILE)
    D = 1.
    L = np.array([50, 100.]) * D
    fdo = None
    ls = np.array([5.]) * D
    dth = np.array([10.]) * PI / 180.
    vf = 0.1
    nc = [2]
    nm = list(range(1, 2+1))
    l_params = assemble_params(L, ls, dth, nc, nm)
    for d_params in l_params:
        d_params['D'] = D 
        d_params['vf'] = vf
        d_params['fdo'] = fdo

    # armo las mallas
    for im, d_params in enumerate(l_params):
        print(d_params)
        mc = Mallacom(**d_params, name=f'malla_{im+1:02d}')
        mc.make_malla()
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
