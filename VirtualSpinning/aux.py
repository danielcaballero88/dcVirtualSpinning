""" modulo de funciones auxiliares
algunas son algebraicas
otras con de manejo de archivos
es una mezcla que resulta util tener por separado """

import numpy as np
from scipy import stats


def find_string_in_file(fid, target, mandatory=False):
    fid.seek(0)  # rewind
    target = target.lower()
    len_target = len(target)
    for linea in fid:
        if target == linea[:len_target].lower():
            ierr = 0
            break
    else:
        if mandatory:
            raise EOFError("final de archivo sin encontrar target string")
        ierr = 1
    return ierr


def iguales(x, y, tol=1.0e-8):
    """ returns True if x is equal to y, both floats """
    return np.abs(x - y) < tol


def calcular_longitud_de_segmento(r0, r1):
    """ dados dos nodos de un segmento
    r0: coordenadas "xy" del nodo 0
    r1: coordenadas "xy" del nodo 1
    calcula el largo del segmento """
    dx = r1[0] - r0[0]
    dy = r1[1] - r0[1]
    return np.sqrt(dx * dx + dy * dy)


def calcular_angulo_de_segmento(r0, r1):
    """ dados dos nodos de un segmento
    r0: coordenadas "xy" del nodo 0
    r1: coordenadas "xy" del nodo 1
    calcula el angulo que forman con el eje horizontal """
    dx = r1[0] - r0[0]
    dy = r1[1] - r0[1]
    if iguales(dx, 0.):  # segmento vertical
        if iguales(dy, 0.):  # segmento nulo
            raise ValueError("Error, segmento de longitud nula!!")
        elif dy > 0.:  # segmento hacia arriba
            return 0.5 * np.pi
        else:  # segmento hacia abajo
            return -0.5 * np.pi
    elif iguales(dy, 0.):  # segmento horizontal
        if dx > 0.:  # hacia derecha
            return 0.
        else:  # hacia izquierda
            return np.pi
    else:  # segmento oblicuo
        if dx < 0.:  # segundo o tercer cuadrante
            return np.pi + np.arctan(dy / dx)
        elif dy > 0.:  # primer cuadrante
            return np.arctan(dy / dx)
        else:  # cuarto cuadrante
            return 2. * np.pi + np.arctan(dy / dx)


def calcular_interseccion_entre_segmentos(r00, r01, r10, r11):
    """ dadas las coordenadas "xy" de los dos nodos para cada segmento
    r00: nodo 0 del segmento 0
    r01: nodo 1 del segmento 0
    r10: nodo 0 del segmento 1
    r11: nodo 1 del segmento 1
    calcula el punto de interseccion
    devuelve una lista [x_in, y_in] si se intersectan
    devuelve None si no se intersectan """
    # ---
    x00, y00 = r00  # coordenadas "xy" del nodo 0 del segmento 0
    x01, y01 = r01  # coordenadas "xy" del nodo 1 del segmento 0
    x10, y10 = r10  # coordenadas "xy" del nodo 0 del segmento 1
    x11, y11 = r11  # coordenadas "xy" del nodo 1 del segmento 1
    # ---
    # chequeo paralelismo
    theta0 = calcular_angulo_de_segmento(r00, r01)
    theta1 = calcular_angulo_de_segmento(r10, r11)
    paralelos = iguales(theta0, theta1, np.pi * 1.0e-8)
    if paralelos:
        return None  # no hay chance de interseccion porque son paralelos los segmentos
    # ---
    # lejania
    # variables auxiliares
    bot0 = np.minimum(y00, y01)
    rig0 = np.maximum(x00, x01)
    top0 = np.maximum(y00, y01)
    lef0 = np.minimum(x00, x01)
    bot1 = np.minimum(y10, y11)
    rig1 = np.maximum(x10, x11)
    top1 = np.maximum(y10, y11)
    lef1 = np.minimum(x10, x11)
    #
    maxBottom = np.maximum(bot0, bot1)  # nodo de la derecha mas a la izquierda entre los dos segmentos
    minRight = np.minimum(rig0, rig1)  # nodo de la derecha mas a la izquierda entre los dos segmentos
    minTop = np.minimum(top0, top1)  # nodo de la derecha mas a la izquierda entre los dos segmentos
    maxLeft = np.maximum(lef0, lef1)  # nodo de la izquierda mas a la derecha entre los dos segmentos
    # ahora chequeo
    lejos = ((maxLeft - minRight) > 1.0e-8) or ((maxBottom - minTop) > 1.0e-8)
    if lejos:
        return None  # no hay chance de interseccion porque estan lejos los segmentos
    # ---
    # interseccion
    # si no son paralelos y si tienen cercania entonces encuentro la interseccion
    # usando un sistema de referencia intrinseco al segmento j0 (chi, eta)
    # calculo los angulos
    theta_rel = theta1 - theta0
    while True:  # ahora me fijo que quede en el rango [0 , 2pi)
        if theta_rel >= 0. and theta_rel < 2. * np.pi:
            break
        elif theta_rel >= 2. * np.pi:
            theta_rel -= 2. * np.pi
        elif theta_rel < 0.:
            theta_rel += 2. * np.pi
    # me fijo que j1 no sea verticales en el sistema intrinseco a j0 (para no manejar pendientes infinitas luego)
    m_rel_inf = iguales(theta_rel, np.pi * 0.5, np.pi * 1.0e-8) or iguales(theta_rel, np.pi * 1.5, np.pi * 1.0e-8)
    # coordenadas en (chi,eta) de los nodos del segmento j1
    C0 = np.cos(theta0)
    S0 = np.sin(theta0)

    def cambio_de_coordenadas(x, y):
        chi = (x - x00) * C0 + (y - y00) * S0
        eta = -(x - x00) * S0 + (y - y00) * C0
        return chi, eta

    chi_01, eta_01 = cambio_de_coordenadas(x01, y01)
    chi_10, eta_10 = cambio_de_coordenadas(x10, y10)
    chi_11, eta_11 = cambio_de_coordenadas(x11, y11)

    # chequeo que el segmento 1 realmente corte al eje del 2
    if np.sign(eta_10) == np.sign(eta_11):
        return None  # si eta no cambia de signo entonces no corta al eje
    # supere varios chequeos
    # calculo valor de chi de la interseccion
    if m_rel_inf:
        chi_in = chi_10
    else:
        m_rel = np.tan(theta_rel)
        chi_in = chi_10 - eta_10 / m_rel
    # ---
    # ahora tengo que estudiar que tipo de interseccion es
    # hay cuatro posibilidades:
    # 1) falsa: por fuera de los segmentos (se cruzan las rectas de los segmentos pero no los segmentos)
    # 2) m-m: los segmentos tienen un punto de interseccion que no se corresponde con ningun extremo de segmento
    # 3) m-e o e-m: la interseccion se da en un extremo de uno de los nodos
    # aqui hay que indicar cual es el nodo de cual segmento [s,n] s es el indice del segmento y n, dentro del segmento, si es el nodo 0 o el 1
    # 4) e-e: la interseccion se da en un punto en el que coinciden un extremo de cada segmento
    # aqui tambien hay que indicar algo, cuales son los nodos [s0,n], [s1,n]
    # en resumen para poder incluir todas las opciones debo devolver varios valores
    # primero las coordenadas del nodo interseccion,
    # segundo el tipo de interseccion (1, 2, 3 o 4),
    # y tercero y cuarto, dos enteros indicando si se corresponde o no con el extremo de los nodos
    # si no es con ninguno de los extremos de algun segmento, en lugar de n=0 o n=1, mando n=None y listo
    # si la interseccion fue por fuera del segmento j0 entonces no es valida
    # ---
    # chequeo posibilidad 1
    dl_0 = chi_01
    fuera_de_seg_0 = (chi_in < 0.) or (chi_in > dl_0)
    if fuera_de_seg_0:
        return None
    # en este punto ya es seguro que hay interseccion
    # calculo x e y:
    x_in = x00 + chi_in * np.cos(theta0)
    y_in = y00 + chi_in * np.sin(theta0)
    r_in = [x_in, y_in]
    # ---
    # posibilidades 3 y 4
    # ahora me fijo si la interseccion coincide con algun extremo
    nin_s0 = None
    nin_s1 = None
    #
    if iguales(chi_in, 0.):
        nin_s0 = 0
    elif iguales(chi_in, dl_0):
        nin_s0 = 1
    #
    if iguales(x_in, x10) and iguales(y_in, y10, tol=1.0e-2):
        nin_s1 = 0
    elif iguales(x_in, x11) and iguales(y_in, y11, tol=1.0e-2):
        nin_s1 = 1
    # ---
    # posibilidad 4 (extremo-extremo)
    if nin_s0 is not None and nin_s1 is not None:  # aqui ninguno es None
        return r_in, 4, nin_s0, nin_s1
    # ---
    # posibilidad 3 (extremo-medio)
    if nin_s0 is None or nin_s1 is None:  # al menos uno es None
        if nin_s0 is not None or nin_s1 is not None:  # al menos uno no es None
            return r_in, 3, nin_s0, nin_s1  # por lo tanto tengo un extremo y uno no (uno es None y otro no lo es)
    # si los dos son None entonces paso a la posibilidad 2
    # ---
    # chequeo posibilidad 2
    # en este punto ya se supone que la interseccion no ha sido ni fuera del intervalo, ni en un extremo
    # igual chequeo por seguridad
    if (chi_in > 0.) and (chi_in < dl_0):
        # la interseccion se da dentro del segmento
        return r_in, 2, None, None
    # ---
    # si llegue hasta aca hay algun error
    raise ValueError("Hay algo que esta muy mal")


def compute_from_curve(x, xx, yy, extrapolar=False):
    """ hay una funcion np.interpolate pero no le confio """
    # primero chequeo si x cae fuera de intervalo
    if x < xx[0] or x > xx[-1]:
        if extrapolar:
            if x < xx[0]:
                return yy[0]
            else:
                return yy[-1]
    # de lo contrario tengo que calcularlo dentro de intervalo interpolando linealmente
    nx = len(xx)
    for i in range(1, nx):
        if x < xx[i]:
            slope = (yy[i] - yy[i - 1]) / (xx[i] - xx[i - 1])
            return yy[i - 1] + slope * (x - xx[i - 1])
