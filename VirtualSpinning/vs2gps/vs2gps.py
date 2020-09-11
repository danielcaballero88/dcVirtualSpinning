"""
Module to read a Mallasim object and write the mesh and param GPsolver files
"""


import numpy as np
from VirtualSpinning.Mallasim.Mallasim import Mallasim 

FMT_INT = "10d"
FMT_FLT = "+17.8e"

def write_mesh(ms, Meshfile): 
    # ---
    # Write Header Data
    Meshfile.write("*NODAL DOFs\n")
    idoft = 4
    Meshfile.write(f"{idoft:{FMT_INT}}\n\n")
    Meshfile.write("*DIMEN\n")
    ndim = 2
    Meshfile.write(f"{ndim:{FMT_INT}}\n\n")
    # ---
    # Write Coordinates
    Meshfile.write("*COORDINATES\n")
    Meshfile.write(f"{ms.nodos.n:{FMT_INT}}\n\n")
    for i in range(ms.nodos.n): 
        r = ms.nodos.r0[i]
        linea = "".join([f"{v:{FMT_FLT}}" for v in r]) + "\n"
        Meshfile.write(linea)
    Meshfile.write("\n")
    # ---
    # Write Groups
    Meshfile.write("*ELEMENT GROUPS\n")
    # header
    ngroups = 1
    Meshfile.write(f"{ngroups:{FMT_INT}}\n\n")
    ifirst = 1
    ilast = ms.fibras.n + ms.nodos.mf.sum()
    Meshfile.write(f"{ifirst:{FMT_INT}}{ilast:{FMT_INT}}{'Generic':>10s}\n\n")
    # groups of fibers
    nodel = 2
    for i in range(ms.fibras.n): 
        Meshfile.write(f"{nodel:{FMT_INT}}\n")
    Meshfile.write("\n")
    # groups of border nodes (dirichlet nodes)
    nodel = 1
    numbornods = ms.nodos.mf.sum()
    for i in range(numbornods): 
        Meshfile.write(f"{nodel:{FMT_INT}}\n")
    Meshfile.write("\n")
    # ---
    # Write Incidence 
    Meshfile.write("*INCIDENCE\n\n")
    # fibras 
    for i in range(ms.fibras.n): 
        n0, n1 = ms.fibras.con[i]
        Meshfile.write("".join([f"{v:{FMT_INT}}" for v in (n0+1,n1+1)]) + "\n")
    Meshfile.write("\n")
    # nodos frontera (dirichlet)
    all_inods = np.arange(ms.nodos.n) + 1
    bor_inods = all_inods[ms.nodos.mf]
    for inod in bor_inods: 
        Meshfile.write(f"{inod:{FMT_INT}}\n")
    Meshfile.write("\n")
    # ---
    # Write Element Type
    Meshfile.write("*ELEMENT TYPE\n\n")
    # fibers
    tipo = 1
    for i in range(ms.fibras.n): 
        Meshfile.write(f"{tipo:{FMT_INT}}\n")
    Meshfile.write("\n")
    # border nodes (dirichlet nodes)
    tipo = 2
    numbornods = ms.nodos.mf.sum()
    for i in range(numbornods): 
        Meshfile.write(f"{tipo:{FMT_INT}}\n")
    Meshfile.write("\n")
    # --- 
    # Write Element Mat
    Meshfile.write("*ELEMENT MAT\n\n")
    # fibers
    for i in range(ms.fibras.n): 
        Meshfile.write(f"{i+1:{FMT_INT}}\n")
    Meshfile.write("\n")
    # border nodes (dirichlet nodes, one group for all)
    numbornods = ms.nodos.mf.sum()
    j = ms.fibras.n + i + 1
    for i in range(numbornods): 
        Meshfile.write(f"{j:{FMT_INT}}\n")
    Meshfile.write("\n")
    # dirichlet conditions (not given here at the moment)
    maskD = np.array([0, 0, 0, 0], dtype=int)
    valsD = np.array([0., 0., 0., 0.], dtype=float)
    # ---
    # Write Dirichlet Conditions
    Meshfile.write("*DIRICHLET CONDITIONS\n\n")
    for i in range(ms.nodos.n * 2): 
        Meshfile.write("".join([f"{v:{FMT_INT}}" for v in maskD]) + "\n")
    Meshfile.write("\n")
    for i in range(ms.nodos.n * 2): 
        Meshfile.write("".join([f"{v:{FMT_FLT}}" for v in valsD]) + "\n")
    Meshfile.write("\n")
    # --- 
    # Done
    Meshfile.write("*END")
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
