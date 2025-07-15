# Collection of functions to perform analysis on H-AMR data

import numpy as np
from .snapshot import Snapshot
import pp_c2 as pp_c

# Calculates various quantities 
def misc_calc(s, calc_bu=1, calc_bsq=1,calc_eu = 0, calc_esq=0):
    if(calc_bu==1):
        bu=np.copy(s.uu)
    else:
        bu=np.zeros((1, 1, 1, 1, 1), dtype=s.rho.dtype)
    if (calc_eu == 1):
        eu = np.copy(s.uu)
    else:
        eu = np.zeros((1, 1, 1, 1, 1), dtype=s.rho.dtype)
    if (calc_bsq == 1):
        bsq=np.copy(s.rho)
    else:
        bsq=np.zeros((1, 1, 1, 1), dtype=s.rho.dtype)
    if (calc_esq == 1):
        esq = np.copy(s.rho)
    else:
        esq = np.zeros((1, 1, 1, 1), dtype=s.rho.dtype)

    pp_c.misc_calc(
        s.bs1new, s.bs2new, s.bs3new, s.nb, s.axisym, s.uu, s.B, s.E, 
        s.bu, s.eu, s.gcov, s.bsq, s.esq, 
        s.calc_bu, s.calc_eu, s.calc_bsq, s.calc_esq)


# Calculate the stress energy tensor T^kappa_nu
def Tcalcud_new(s, kappa, nu):
    bd_nu = (s.gcov[nu,:]*s.bu).sum(0)
    ud_nu = (s.gcov[nu,:]*s.uu).sum(0)
    Tud = s.bsq * s.uu[kappa] * ud_nu + 0.5 * s.bsq * (kappa==nu) - s.bu[kappa] * bd_nu +(s.rho + s.ug + (s.gam - 1) * s.ug) * s.uu[kappa] * ud_nu + (gam - 1) * ug * (kappa==nu)
    return Tud

# Calculate the total mass of disk in code units
def calc_Mtot(s):
    return np.sum((s.rho * s.uu[0]) * s._dx1 * s._dx2 * s._dx3 * s.gdet)

# Calculate mass accretion rate as function of radius
def calc_Mdot(s):
    return (-s.gdet * s.rho * s.uu[1] * s._dx2 * s._dx3).sum(-1).sum(-1)

# Calculate energy accretion rate as function of radius
def calc_Edot(s):
    temp = Tcalcud_new(s, 1, 0)* s.gdet * s._dx2 * s._dx3
    Edot = (temp).sum(-1).sum(-1)
    Edotj = (temp*(s.bsq/s.rho>3)).sum(-1).sum(-1)
    return Edot, Edotj

def calc_Ldot(s):
    return (Tcalcud_new(s, 1, 3)* s.gdet * s._dx2 * s._dx3).sum(-1).sum(-1)

# Calculate magnetic flux phibh as function of radius
def calc_phibh(s):
    phibh = 0.5 * (np.abs(s.gdet * s.B[1]) * s._dx2 * s._dx3).sum(-1).sum(-1)

def calc_rad_avg(s):
    return (s.r * s.rho * s.gdet * s._dx1 * s._dx2 * s._dx3).sum(-1).sum(-1).sum(-1) / (
        (s.rho * s.gdet * s._dx1 * s._dx2 * s._dx3).sum(-1).sum(-1).sum(-1))


