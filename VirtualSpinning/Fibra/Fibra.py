import math
import numpy as np
from VirtualSpinning.aux import find_string_in_file


class Fibra(object):
    """
    Clase para calcular la curva de una sola fibra traccionada
    """
    def __init__(self, Et, EbEt, doteps, s0, nh, lamr=1., tenbrk=1000.):
        self.param = {
            'Et' : Et,
            'Eb' : EbEt * Et,
            'doteps' : doteps,
            's0' : s0,
            'nh' : nh,
            'tenbrk': tenbrk
        }
        self.broken = False
        self.lamr = lamr
        self.lamp = 1.

    @classmethod 
    def from_cf(cls, cf, lamr=1.): 
        with open(cf, 'r') as f:
            find_string_in_file(f, "* Parametros constitutivos") 
            _ = int(next(f)) 
            param = [float(val) for val in next(f).replace('d','e').split()]
        ilaw, Et, EbEt, doteps, s0, nh, tenbrk = param
        assert(ilaw == 4) 
        self = cls(Et, EbEt, doteps, s0, nh, lamr, tenbrk)
        return self

    def calc_ten(self, lam):
        """ 
        Calcula la tension en un incremento de tiempo
        puede haber plasticidad y puede estar rota
        pero aca no se incrementan esas variables
        """
        # Tomo algunas variables mas comodas
        Et = self.param['Et']
        Eb = self.param['Eb']
        lamrp = self.lamr * self.lamp
        # Calculo segun el caso
        if self.broken:  # fibra rota
            ten = 0.
        elif lam < lamrp:  # fibra enrulada
            ten = Eb * (lam - 1.)
        else:  # fibra reclutada
            tenr = Eb*(lamrp - 1.)  # tension en el punto de reclutamiento
            ten = tenr + Et * (lam / lamrp - 1.)
        return ten

    def calc_plas(self, ten, dt):
        """
        Calcula la tasa de deformacion plastica en funcion de la tension
        Tambien la rotura si se produce
        """

        # Calculo la tasa de plasticidad y/o si rompe la fibra
        if self.broken:
            dotlamp = 0.
        elif ten > self.param['tenbrk']:
            self.broken = True 
            dotlamp = 0. 
        else: 
            s = self.param['s0'] * self.lamp**self.param['nh']
            dotlamp = self.param['doteps'] * math.sinh(ten / s)

        # Incremento la plasticidad
        self.lamp = self.lamp + dotlamp * dt

    def traccionar(self, dt, dotlam, lamf):
        """
        Calcular la curva tension vs lam para una traccion en el tiempo
        """
        # Variables iniciales del esquema temporal
        time = 0.
        lam = 1.
        # Listas de variables principales a guardar
        rec_time = [time]
        rec_lam = [lam]
        rec_ten = [0.]
        # Lista de otras variables a guardar
        rec_lamp = [self.lamp] 
        rec_lam_ef = [1. / self.lamr / self.lamp]
        while lam < lamf:
            time += dt 
            lam += dotlam * dt
            ten = self.calc_ten(lam)
            self.calc_plas(ten, dt)
            rec_time.append(time) 
            rec_lam.append(lam) 
            rec_ten.append(ten)
            rec_lamp.append(self.lamp) 
            rec_lam_ef.append(lam / self.lamr / self.lamp)
        rec = {
            'time': np.array(rec_time),
            'lam': np.array(rec_lam),
            'eps': np.array(rec_lam) - 1.,
            'ten': np.array(rec_ten),
            'lamp': np.array(rec_lamp), 
            'epsp': np.array(rec_lamp) - 1.,
            'lam_ef': np.array(rec_lam_ef),
            'eps_ef': np.array(rec_lam_ef) -1.
        }
        return rec
