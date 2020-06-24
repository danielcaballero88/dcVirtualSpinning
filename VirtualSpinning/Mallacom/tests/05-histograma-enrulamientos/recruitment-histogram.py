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
    mallas = [] 
    for i in range(4):
        mc = Mallacom.leer_de_archivo(DIR / f'malla_{i+1:02d}.txt')
        mallas.append(mc)

    # Grafico los histogramas por separad
    for i in range(4):
        fig, ax = plt.subplots()
        mallas[i].graficar_histograma_enrulamientos(fig, ax)

    # Grafico todos en un solo grafico
    fig, ax = plt.subplots()
    colores = ['r', 'g', 'b', 'y']
    for i in range(4):
        mallas[i].graficar_histograma_enrulamientos(
            fig, ax, dx_mult=0.25, x_offset=-0.375+i*0.25, 
            lamr_min=1.0, lamr_max=1.1, nbins=10, 
            color=colores[i], label=str(i+1)
        )
    ax.legend() 

    plt.show()


if __name__ == '__main__':
    main()
