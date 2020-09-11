from VirtualSpinning.Mallasim.Mallasim import Mallasim 
from VirtualSpinning.vs2gps import vs2gps
import os 

cwd = os.path.dirname(__file__)
mspath = os.path.join(cwd, "malla_i_s.txt")
ms = Mallasim.leer_de_archivo(mspath)

if False:
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots()
    ms.pre_graficar_bordes(fig, ax)
    colores = ['blue', 'red']
    ms.pre_graficar(fig, ax, cby="lamr",
                    cbar=True, cmap="rainbow",
                    plot_afin=False)
    figpath = os.path.join(cwd, "Mallasim.png")
    fig.savefig(figpath)


meshpath = os.path.join(cwd, "Mesh.txt")
with open(file=meshpath, mode="w") as f: 
    vs2gps.write_mesh(ms, f)

parampath = os.path.join(cwd, "Param.txt")
with open(file=parampath, mode="w") as f:
    vs2gps.write_param(ms,f)