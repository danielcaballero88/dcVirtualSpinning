from numpy import pi as PI
from matplotlib import pyplot as plt
from VirtualSpinning.Mallacom import Mallacom


def main():
    print('Hi!, I\'m test basic!')
    Dm = 1.
    L = 50. * Dm
    fundisor = None
    dl = 5. * Dm
    devang = 10. * PI / 180.
    volfrac = 0.1
    ncaps = 2
    mc = Mallacom(L, Dm, volfrac, dl, devang, fundisor=None)
    for i in range(1, ncaps + 1):
        mc.make_capa(dl, Dm, devang, volfrac, orient_distr=fundisor)
    mc.guardar_en_archivo('malla_test_basic.txt')
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    mc.pre_graficar_bordes(fig, ax)
    mc.pre_graficar_fibras(fig, ax, color_por="nada", byn=True, linewidth=1.5)
    plt.show()

if __name__ == '__main__':
    main()
