from VirtualSpinning.Mallasim.Mallasim import Mallasim 
from VirtualSpinning.vs2gps import vs2gps
import os 

cwd = os.path.dirname(__file__)
mspath = os.path.join(cwd, "malla_i_s.txt")
ms = Mallasim.leer_de_archivo(mspath)

meshpath = os.path.join(cwd, "Mesh.txt")
with open(file=meshpath, mode="w") as f: 
    vs2gps.write_mesh(ms, f)

parampath = os.path.join(cwd, "Param.txt")
with open(file=parampath, mode="w") as f:
    vs2gps.write_param(ms,f)