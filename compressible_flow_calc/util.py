# General aerodynamic utility functions
# Distrubuted as a part of compressible_flow_calc, under GNU GPL3.
# (C) Ariel Mordoch 2018
import numpy as np
from . import unit, Q


def v(M, a):
    """Velocity, given Mach num."""
    return M * a


def M(v, a):
    """Mach num., given speed of sound and velocity."""
    return v / a


def a(gamma, R, T):
    """Speed of sound, given gamma, R_specific, and temperature."""
    return np.sqrt(gamma * R * T)


def Re(V, L, rho=Q(1.225, 'kg/m^3'), mu=Q(1.789e-5, "Pa*s"), nu=None):
    """Reynolds Number, given rho, V, and L, as well as either kinematic or
    dynamic viscosity. If only V and L are given, the Re at sea level standard
    conditions is given."""
    if nu is None:
        return (rho * V * L) / mu
    else:
        return (V * L) / nu

