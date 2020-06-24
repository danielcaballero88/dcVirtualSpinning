"""
En este test leo de archivos y grafico histogramas de orientaciones
de dos Mallacom, una con fdo uniforme y una mas alineada en pi/2
"""


from matplotlib import pyplot as plt
from pathlib import Path
from VirtualSpinning.Mallacom.Mallacom import Mallacom


CWD = Path.cwd()
FILE = Path(__file__)
DIR = FILE.parent


def main():
    print('testing: Mallacom -> read-plot')
    print('CWD: ', CWD)
    print('DIR: ', DIR)
    print('__file__: ', FILE)

    # Leo las mallas
    mc_u = Mallacom.leer_de_archivo(DIR / 'malla_uniforme.txt')
    mc_a = Mallacom.leer_de_archivo(DIR / 'malla_alineada.txt')

    # Grafico el histograma de la malla uniforme
    fig, ax = plt.subplots()
    mc_u.graficar_histograma_orientaciones(fig, ax)

    # Grafico el histograma de la malla alineada
    fig, ax = plt.subplots()
    mc_a.graficar_histograma_orientaciones(fig, ax)

    # Grafico los dos en uno solo 
    fig, ax = plt.subplots()
    mc_u.graficar_histograma_orientaciones(fig, ax, dx_mult=0.5, x_offset = -0.25, label='u') 
    mc_a.graficar_histograma_orientaciones(fig, ax, dx_mult=0.5, x_offset = +0.25, color='white', label='a')
    ax.legend() 

    plt.show()


if __name__ == '__main__':
    main()
