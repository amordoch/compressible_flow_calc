# Compressible flow calculator - iteration functions
# By Ariel Mordoch - distributed under GNU GPL3

import math
from numpy import linspace
from .calc import *

# Tuples of iteration functions
_ISENFUNCS = (A_over_Astar, T0_over_T, p0_over_p,
              rho0_over_rho)
_SHOCKFUNCS = (p02_over_p01, T2_over_T1, p2_over_p1)
_EXPANSIONFUNCS = (nu)


def iterate(correct_output, x, y, accuracy):
    """
    iterate() takes in a given function output [y] = [f(x)],
    a list of input values [x], and a list of output values [f(x)],
    and iterates through the list to find
    where the desired output occurs.

    :param float correct_output: known correct function output (y value)
    :param list x: list of input values
    :param list y: list of output values
    :param float accuracy: desired percent error
    :return: x value found, or None for no solution
    """
    accuracy = abs(accuracy)
    i = 0
    x_val = None
    for output in y:
        pct_err = abs((correct_output - output) / correct_output)
        if pct_err <= accuracy:
            x_val = x[i]
            break
        else:
            i += 1
    return x_val


# Iteration over isentropic functions #
def isentropic_iterate(gamma, func, correct_output, start_M=1, end_M=5, accuracy=1e-4, step=1e-5):
    """ Iterates over isentropic flow functions using input step size and accuracy """
    # Create mach num list
    # + 1 to include end_M as the final step
    mach_nums = list(linspace(start_M, end_M, (end_M - start_M) / step + 1))
    func = _ISENFUNCS[func]
    # Generate list of the function's outputs
    # Code will break if all functions do not take uniform inputs
    outputs = [func(M, gamma) for M in mach_nums]
    # Finally, iterate.
    found_M = iterate(correct_output, mach_nums, outputs, accuracy)
    if found_M is not None:
        # Round found_M to the specified significance
        return round(found_M, abs(int(math.log10(accuracy))))
    else:
        return None


# Iteration over normal shock functions #
def shock_iterate(gamma, func, correct_output, startM=1, endM=5, accuracy=1e-4, step=1e-5):
    """ Iterates over normal shock functions using input step size and accuracy. Modeled after isentropic_iterate() """
    # Create mach num list
    # + 1 to include end_M as the final step
    mach_nums = list(linspace(startM, endM, (endM - startM) / step + 1))
    func = _SHOCKFUNCS[func]
    # Generate list of outputs of given function and iterate
    outputs = [func(M1, M2(M1, gamma), gamma) for M1 in mach_nums]
    found_M = iterate(correct_output, mach_nums, outputs, accuracy)
    if found_M is None:
        return None
    else:
        # Return iterated value terminated to the desired accuracy
        return round(found_M, abs(int(math.log10(accuracy))))


# Iteration over expansion functions #
def expansion_iterate(gamma, func, correct_output, startM=1, endM=5, accuracy=1e-4, step=1e-5):
    """ Iterates over nu(M, gamma) to find the mach number """
    # Create mach num list
    # + 1 to include end_M as the final step
    mach_nums = list(linspace(startM, endM, (endM - startM) / step + 1))
    func = _EXPANSIONFUNCS  # just calc.nu() for now
    # Generate list of nu values and iterate
    nu_vals = [func(M, gamma) for M in mach_nums]
    found_M = iterate(correct_output, mach_nums, nu_vals, accuracy)
    if found_M is None:
        return None
    else:
        # Return value with the given accuracy
        return round(found_M, abs(int(math.log10(accuracy))))
