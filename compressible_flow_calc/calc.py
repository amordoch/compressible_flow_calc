# Compressible flow calculator - utility functions
# By Ariel Mordoch - distributed under GNU GPL3 
import math


# Isentropic flow
def A_over_Astar(M, gamma):
    return 1 / M * ((1 + (gamma - 1) / 2 * M ** 2) / ((gamma + 1) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))


def mDot_over_A(M, gamma, p0, R, T0):
    const = p0 / math.sqrt(R * T0)
    f = math.sqrt(gamma) * M / (1 + (gamma - 1) / 2 * M ** 2) ** ((gamma + 1) / (2 * (gamma - 1)))
    return const * f


def T0_over_T(M, gamma):
    return 1 + (gamma - 1) / 2 * M ** 2


def p0_over_p(M, gamma):
    return (T0_over_T(M, gamma)) ** (gamma / (gamma - 1))


def rho0_over_rho(M, gamma):
    return (T0_over_T(M, gamma)) ** (1 / (gamma - 1))


# Normal shocks
def mu(M):
    return math.asin(1 / M)


def M2(M1, gamma):
    numer = M1 ** 2 + 2 / (gamma - 1)
    denom = 2 * gamma / (gamma - 1) * M1 ** 2 - 1
    return math.sqrt(numer / denom)


def T2_over_T1(M1, M2, gamma):
    T02_T2 = T0_over_T(M2, gamma)
    T01_T1 = T0_over_T(M1, gamma)
    return (T02_T2) ** -1 / (T01_T1) ** -1


def v2_over_v1(M1, M2, gamma):
    return M2 / M1 * math.sqrt(T2_over_T1(M1, M2, gamma))


def rho2_over_rho1(M1, M2, gamma):
    return (v2_over_v1(M1, M2, gamma)) ** -1


def p2_over_p1(M1, M2, gamma):
    return (1 + gamma * M1 ** 2) / (1 + gamma * M2 ** 2)


def p02_over_p01(M1, M2, gamma):
    p02_p2 = p0_over_p(M2, gamma)
    p01_p1 = p0_over_p(M1, gamma)
    return p02_p2 / p01_p1 * p2_over_p1(M1, M2, gamma)


def p02_over_p1(M1, M2, gamma):
    return p0_over_p(M2, gamma) * p2_over_p1(M1, M2, gamma)


def rho02_over_rho01(M1, M2, gamma):
    return p02_over_p01(M1, M2, gamma)


def A2star_over_A1star(M1, M2, gamma):
    return (p02_over_p01(M1, M2, gamma) ** -1)


# Oblique shocks

def beta(M1, gamma, theta, alpha=1):
    # Compute this elementary function
    L1 = (M1 ** 2 - 1) ** 2
    L2 = 3 * (1 + (gamma - 1) / 2 * (M1 ** 2)) * (1 + (gamma + 1) / 2 * (M1 ** 2)) * (math.tan(theta)) ** 2
    Lambda = math.sqrt(L1 - L2)
    # Cube it, then compute whatever this is
    chi = 1 / Lambda ** 3 * ((M1 ** 2 - 1) ** 3 - 9 * (1 + (gamma - 1) / 2 * M1 ** 2) * (
                1 + (gamma - 1) / 2 * M1 ** 2 + (gamma + 1) / 4 * M1 ** 4) * (math.tan(theta)) ** 2)
    # Good god...
    tanB = (M1 ** 2 - 1 + 2 * Lambda * math.cos((4 * math.pi * alpha + math.acos(chi)) / 3)) / (
                3 * (1 + (gamma - 1) / 2 * M1 ** 2) * math.tan(theta))
    # Someone much smarter than me found this formula
    beta = math.atan(tanB)
    return beta


def M1n(M1, beta):
    return M1 * math.sin(beta)


def M2n_M2(M2n, beta, theta):
    return M2n / math.sin(beta - theta)


# Prandtl-Meyer expansions
def nu(M, gamma):
    return math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(
        math.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) - math.atan(math.sqrt(M ** 2 - 1))


def normal_shock_final(M1, gamma):
    """
    Solves for flow property ratios after a normal shock (in no particular
    order)
    :param float M1: mach number before shock
    :param float gamma: ratio of specific heats (C_p / C_v)
    :return: a dict of property ratios, with division replaced by underscores:
    i.e. 'T2 / T1 -> T2_T1'
    :rtype: dict
    """

    M_final = M2(M1, gamma)
    return {
        'M2': M_final,
        'T2_T1': T2_over_T1(M1, M_final, gamma),
        'rho2_rho1' : rho2_over_rho1(M1, M_final, gamma),
        'p2_p1': p2_over_p1(M1, M_final, gamma),
        'p02_p01':  p02_over_p01(M1, M_final, gamma),
        'p02_p1': p02_over_p1(M1, M_final, gamma),
        'rho02_rho01': rho02_over_rho01(M1, M_final, gamma),
        'A2star_A1star': A2star_over_A1star(M1, M_final, gamma),
        'v2_v1': v2_over_v1(M1, M_final, gamma)
    }


def oblique_shock_final(M1, theta, gamma, alpha=1):
    """
    Solves for flow property ratios after an oblique shock (in no particular
    order)
    :param float M1: mach number before shock
    :param float theta: turning angle
    :param float gamma: ratio of specific heats (C_p / C_v)
    :param int alpha: 0 for a strong oblique shock, 1 for a weak oblique shock
    :return: a dict of property ratios, same format as normal_shcck_format()
    :rtype: dict
    """

    shock_angle = beta(M1, gamma, theta, alpha)
    M1_normal = M1n(M1, shock_angle)
    props = normal_shock_final(M1_normal, gamma)
    # M2 is actually going to be M2n, must correct this
    props['M2n'] = props['M2']
    props['M2'] = M2n_M2(props['M2'], shock_angle, theta)
    # For oblique shocks, p02/p1 is invalid, so we can remove it
    props.pop('p02_p1')
    # Add M1n and beta to output
    props['M1n'] = M1_normal
    props['beta'] = shock_angle
    return props

