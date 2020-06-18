from matplotlib import pyplot as plt
from pathlib import Path
from VirtualSpinning.Mallasim.Mallasim import Mallasim
from VirtualSpinning.Mallacom.Mallacom import Mallacom


CWD = Path.cwd()
FILE = Path(__file__)
DIR = FILE.parent


def main():
    print('testing: Mallacom -> read-plot')
    print('CWD: ', CWD)
    print('DIR: ', DIR)
    print('__file__: ', FILE)
    mc = Mallacom.leer_de_archivo(DIR / 'malla_test_basic.txt')
    ms = Mallasim.leer_desde_archivo(DIR / 'malla_test_basic_i_s.txt')
    fig, ax = plt.subplots()
    mc.pre_graficar_fibras(fig, ax)
    ms.pre_graficar_bordes(fig, ax)
    colores_cm = ['blue', 'red']
    ms.pre_graficar(fig, ax, color_por="lamr", linewidth=1.5,
                    # lam_min=0.0, lam_max=100.,
                    barracolor=True, colormap="rainbow", colores_cm=colores_cm, maxnfibs=3000,
                    afin=False, colorafin="k", linewidthafin=1.5)
    plt.show()


if __name__ == '__main__':
    main()
