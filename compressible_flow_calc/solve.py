import pint
import numpy as np
from compressible_flow_calc import unit, Q_

units = pint.UnitRegistry()
# Some constants
R_air = 287 * units.joule / units.mol / units.K

def oblique_solve(gamma, **kwargs):
    pass


def normal_solve(gamma, **kwargs):
    # Determine what properties we have and use them to get initial mach number
    if 'M1' in kwargs.keys():
        M1 = kwargs['M1']
    elif 'v1' in kwargs and 'T1' in kwargs and 'R' in kwargs:
        # Must use numpy as math is not compatible with pint
        a = np.sqrt(gamma * kwargs['R'] * kwargs['T1'])
        M1 = kwargs['v1'] * a

