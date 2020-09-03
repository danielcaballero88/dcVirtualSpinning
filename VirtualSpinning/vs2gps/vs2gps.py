"""
Module to read a Mallasim object and write the mesh and param GPsolver files
"""


import numpy as np
from VirtualSpinning.Mallasim.Mallasim import Mallasim 

FMT_INT = "10d"
FMT_FLT = "+17.8e"

def write_mesh(ms, mf): 
    # ---
    # Write Header Data
    mf.write("*NODAL DOFs\n")
    idoft = 2
    mf.write(f"{idoft:{FMT_INT}}\n\n")
    mf.write("*DIMEN\n")
    ndim = 2
    mf.write(f"{ndim:{FMT_INT}}\n\n")
    # ---
    # Write Coordinates
    mf.write("*COORDINATES\n")
    mf.write(f"{ms.nodos.n:{FMT_INT}}\n\n")
    for i in range(ms.nodos.n): 
        r = ms.nodos.r0[i]
        linea = "".join([f"{v:{FMT_FLT}}" for v in r]) + "\n"
        mf.write(linea)
    mf.write("\n")
    # ---
    # Write Groups
    mf.write("*ELEMENT GROUPS\n")
    # header
    ngroups = 1
    mf.write(f"{ngroups:{FMT_INT}}\n\n")
    ifirst = 1
    ilast = ms.fibras.n + ms.nodos.mf.sum()
    mf.write(f"{ifirst:{FMT_INT}}{ilast:{FMT_INT}}{'Generic':>10s}\n\n")
    # groups of fibers
    nodel = 2
    for i in range(ms.fibras.n): 
        mf.write(f"{nodel:{FMT_INT}}\n")
    mf.write("\n")
    # groups of border nodes (dirichlet nodes)
    nodel = 1
    numbornods = ms.nodos.mf.sum()
    for i in range(numbornods): 
        mf.write(f"{nodel:{FMT_INT}}\n")
    mf.write("\n")
    # ---
    # Write Incidence 
    mf.write("*INCIDENCE\n\n")
    # fibras 
    for i in range(ms.fibras.n): 
        n0, n1 = ms.fibras.con[i]
        mf.write("".join([f"{v:{FMT_INT}}" for v in (n0+1,n1+1)]) + "\n")
    mf.write("\n")
    # nodos frontera (dirichlet)
    all_inods = np.arange(ms.nodos.n) + 1
    bor_inods = all_inods[ms.nodos.mf]
    for inod in bor_inods: 
        mf.write(f"{inod:{FMT_INT}}\n")
    mf.write("\n")
    # ---
    # Write Element Type
    mf.write("*ELEMENT TYPE\n\n")
    # fibers
    tipo = 1
    for i in range(ms.fibras.n): 
        mf.write(f"{tipo:{FMT_INT}}\n")
    mf.write("\n")
    # border nodes (dirichlet nodes)
    tipo = 2
    numbornods = ms.nodos.mf.sum()
    for i in range(numbornods): 
        mf.write(f"{tipo:{FMT_INT}}\n")
    mf.write("\n")
    # --- 
    # Write Element Mat
    mf.write("*ELEMENT MAT\n\n")
    # fibers
    for i in range(ms.fibras.n): 
        mf.write(f"{i+1:{FMT_INT}}\n")
    mf.write("\n")
    # border nodes (dirichlet nodes)
    numbornods = ms.nodos.mf.sum()
    for i in range(numbornods): 
        j = ms.fibras.n + i + 1
        mf.write(f"{j:{FMT_INT}}\n")
    mf.write("\n")
    # dirichlet conditions (not given here at the moment)
    mD = np.array([0, 0], dtype=int)
    uD = np.array([0., 0.], dtype=float)
    # ---
    # Write Dirichlet Conditions
    mf.write("*DIRICHLET CONDITIONS\n\n")
    for i in range(ms.nodos.n): 
        mf.write("".join([f"{v:{FMT_INT}}" for v in mD]) + "\n")
    mf.write("\n")
    for i in range(ms.nodos.n): 
        mf.write("".join([f"{v:{FMT_FLT}}" for v in uD]) + "\n")
    mf.write("\n")
    # --- 
    # Done
    mf.write("*END")
    return 0
    # ---


def write_param(ms, mf):
    # ---
    # Write Parameter Groups
    mf.write("*Parameter Groups\n")
    ngroups = ms.fibras.n + 1
    mf.write(f"{ngroups:{FMT_INT}}\n\n")
    # ---
    # Write Real Parameters
    mf.write("*Real Parameters\n\n")
    # number of parameters per element
    nparamfibs = 2
    for i in range(ms.fibras.n):
        mf.write(f"{nparamfibs:{FMT_INT}}\n")
    mf.write("\n")
    nparamdirich = 1
    mf.write(f"{nparamdirich:{FMT_INT}}\n\n")
    # values per element 
    for i in range(ms.fibras.n):
        lamr = ms.fibras.lamsr[i] 
        d = ms.fibras.ds[i]
        fibparam = [lamr, d]
        mf.write("".join([f"{v:{FMT_FLT}}" for v in fibparam]) + "\n")
    mf.write("\n")
    dirichparam = [0]
    mf.write("".join([f"{v:{FMT_FLT}}" for v in dirichparam]) + "\n\n")
    # --- 
    # Write Integer Parameters 
    mf.write("*Integer Parameters\n\n")
    # just a lot of zeros to dictate no integer parameters is fine 
    nparamfibs = 0
    for i in range(ms.fibras.n):
        mf.write(f"{nparamfibs:{FMT_INT}}\n")
    mf.write("\n")
    nparamdirich = 0
    mf.write(f"{nparamdirich:{FMT_INT}}\n\n")
    # --- 
    # Done
    mf.write("*END")
    return 0
    # ---
