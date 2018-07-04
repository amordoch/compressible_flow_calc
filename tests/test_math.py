# Tests for mathematical functions in calc.py
# Correct values taken from:
# Anderson, John D. Fundamentals of Aerodynamics. 5th ed., McGraw-Hill, 2011.
from compressible_flow_calc import calc

ACCURACY = 1e-3  # limitation of source text
# DO NOT MODIFY OR ALL TESTS WILL FAIL
GAMMA = 1.4

def is_accurate(true, calculated):
    """
    Checks whether calculated is within ACCURACY (default 1e-4) of true.
    :param float true: known accurate value
    :param float calculated: value to test
    :return: whether calculated is within ACCURACY or true
    :rtype: bool
    """
    pct_diff = abs((true - calculated)/true)
    if pct_diff <= ACCURACY:
        return True
    else:
        return False

def test_isentropic():

    # Source : Anderson Appendix A
    # p0/p, rho0/rho, T0/T, A/A*
    # Embedded list values respective to mach_nums
    mach_nums = [.2e-1, .6, 2.4]
    known_good = [[.1e1, .1e1, .1e1, .2894e2],
                  [.1276e1, .1190e1, .1072e1, .1188e1],
                  [.1462e2, .6794e1, .2152e1, .2403e1]]

    check_vals = []
    check_results = []
    z = 0
    for M in mach_nums:
        kn_gd = known_good[z]
        p0_p = calc.p0_over_p(M, GAMMA)
        rho0_0 = calc.rho0_over_rho(M, GAMMA)
        T0_T = calc.T0_over_T(M, GAMMA)
        A_Astar = calc.A_over_Astar(M, GAMMA)
        i = 0
        vals = [p0_p, rho0_0, T0_T, A_Astar]
        for val in vals:
            assert is_accurate(kn_gd[i], val),\
                'interrogated %.8f, correct is %.4f. mach = %.4f, func = %d' \
                % (val, kn_gd[i], M, vals.index(val))
            i += 1
        z += 1
    print('Isentropic functions are accurate to within'
          + ' %f percent of known good values' % (ACCURACY*100))

test_isentropic()