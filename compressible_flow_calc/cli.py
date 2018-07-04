# CLI for compressible flow calculations
# By Ariel Mordoch, distributed under GNU GPL
import math
from compressible_flow_calc import iteration, calc


def prompt(text='?>', expected=str, default=None):
    """ Prompts user with given text with type-checking """
    while True:
        ans = input(text)
        if ans:
            if expected is not str:
                try:
                    ans = expected(ans)
                    break
                except ValueError:
                    print('Error: %s is not %s' % (ans, str(expected)))
                    continue
            else:
                break
        elif not ans:
            # Fall back to default value is no input is given
            if default is None:
                continue
            else:
                ans = default
                break
    return ans


def confirm(text):
    """ Simple confirmation using prompt() """
    if not text:
        text = '[y/N]'
    while True:
        ans = prompt(text, default='n')
        if ans.lower() == 'y':
            ans = True
            break
        elif ans.lower() == 'n':
            ans = False
            break
    return ans


class CompressibleCLI:

    def __init__(self):
        self.__gamma = 1.4
        self.start()

    def start(self):
        # Intro
        print("Compressible Flow Calculator 1.0")
        print("By Ariel Mordoch")
        self.__gamma = prompt("Enter a value for gamma (leave blank for 1.4): ", float, 1.4)
        # Start Looping
        loop = True
        while loop:
            # Decide between iteration & calculations
            print("---\nSelect option: ")
            print("1: Known Mach number")
            print("2: Unknown Mach number")
            print("3: Exit")
            sel = prompt(expected=int)
            if sel == 1:
                self.knownM()
            elif sel == 2:
                self.unknownM()
            else:
                print("Goodbye!")
                loop = False

    def print_row(self, table, M):
        if table == 1:
            # Isentropic print
            A_As = calc.AoverAstar(M, self.__gamma)
            T0_T = calc.T0_over_T(M, self.__gamma)
            p0_p = calc.p0_over_p(M, self.__gamma)
            rho0_rho = calc.rho0_over_rho(M, self.__gamma)
            print("A/A* = %f \nT0/T = %f \np0/p = %f\nrho0/rho = %f" % (A_As, T0_T, p0_p, rho0_rho))
        elif table == 2 or table == 3:
            # Normal shocks/oblique shocks
            if table == 3:
                theta = prompt("Enter turning angle (deg): ", expected=float) * math.pi / 180
                alpha = prompt("For a strong oblique shock, enter 0: ", default=1)
                if alpha == "0":
                    beta = calc.beta(M, self.__gamma, theta, 0)
                else:
                    beta = calc.beta(M, self.__gamma, theta)
                M1_real = M
                M = M * math.sin(beta) # M1n
                M2 = calc.M2(M, self.__gamma)  # M2n
                M2_real = M2 / math.sin(beta - theta)  # M2
            else:
                M2 = calc.M2(M, self.__gamma)

            T2_T1 = calc.T2overT1(M, M2, self.__gamma)
            rho2_rho1 = calc.rho2_over_rho1(M, M2, self.__gamma)
            p2_p1 = calc.p2_over_p1(M, M2, self.__gamma)
            p02_p01 = calc.p02_over_p01(M, M2, self.__gamma)
            p02_p1 = calc.p02_over_p1(M, M2, self.__gamma)
            rho02_rho01 = calc.rho02_over_rho01(M, M2, self.__gamma)
            A2star_A1star = calc.A2star_over_A1star(M, M2, self.__gamma)
            v2_v1 = calc.v2overv1(M, M2, self.__gamma)
            if table == 3:
                # Some things must be modified for oblique shocks
                print("M1 = %f\nM1n = %f\nM2n = %f\nM2 = %f" % (M1_real, M, M2, M2_real))
                print("T2 / T1 = %f\nrho2 / rho1 = %f\np2 / p1 = %f\np02 / p01 = %f" % (T2_T1, rho2_rho1, p2_p1,
                                                                                        p02_p01))
                print("rho02 / rho01 = %f\nA2* / A1* = %f\nv2 / v1 = %f" % (rho02_rho01, A2star_A1star, v2_v1))
            else:
                print("M2 = %f\nT2 / T1 = %f\nrho2 / rho1 = %f\np2 / p1 = %f\np02 / p01 = %f" % (
                    M2, T2_T1, rho2_rho1, p2_p1, p02_p01))
                print("p02 / p1 = %f\nrho02 / rho01 = %f\nA2* / A1* = %f\nv2 / v1 = %f" % (
                    p02_p1, rho02_rho01, A2star_A1star, v2_v1))
        elif table == 4:
            # Expansions
            nu = calc.nu(M, self.__gamma) * 180 / math.pi
            mu = calc.mu(M) * 180 / math.pi
            print("nu = %f deg\nmu = %f deg" % (nu, mu))

    def knownM(self):
        # Select isentropic, normal shock, oblique shock, or expansion
        print("Select option:")
        print("1: Isentropic Flow Properties\n2: Normal Shock Properties")
        print("3: Oblique Shock Properties\n4: Prandtl-Meyer\n5: Back")
        # Inputs
        sel = prompt(expected=int)
        # mach number
        if sel < 5:
            M = prompt(text="Enter Mach number: ", expected=float)
            # error handling required for mathematical exceptions
            try:
                self.print_row(sel, M)
            except Exception as e:
                print("Error: ", e)
                print("Continuing...")
        else:
            pass

    def unknownM(self):
        print('---\nIterate over:')
        print('1: Isentropic flow values\n2: Normal shock values\n3: Prandtl-Meyer values')
        sel = prompt(expected=int)
        if sel == 1:
            self.isentropic_iter_wrapper()
        elif sel == 2:
            self.shock_iter_wrapper()
        elif sel == 3:
            self.expansion_iter_wrapper()
        else:
            pass

    def isentropic_iter_wrapper(self):
        """ Command line wrapper over iteration.isentropic_iterate()"""
        print("Select function to iterate over:")
        print("1: A / A*")
        print("2: T0 / T")
        print("3: p0 / p")
        print("4: rho0 / rho")
        func = prompt(expected=int)
        settings = self.mbound_err_prompt()
        correct_output = prompt('Enter the function\'s correct output: ', expected=float)
        print("Working...calculation may take some time depending on size of interval and desired accuracy")
        if settings is not None:
            # Correct func choice for 0 indexing
            result = iteration.isentropic_iterate(self.__gamma, func - 1, correct_output, settings[0], settings[1], settings[2],
                                                  settings[3])
        else:
            result = iteration.isentropic_iterate(self.__gamma, func - 1, correct_output)
        if result is not None:
            print('M = %f' % result)
        else:
            print("No solution found within the specific interval")

    def shock_iter_wrapper(self):
        """ Command line wrapper over iteration.shockIteratre() """
        print("Select function to iterate over:")
        print("1: p02 / p01")
        print("2: T2 / T1")
        print("3: p2 / p1")
        func = prompt(expected=int)
        settings = self.mbound_err_prompt()
        correct_output = prompt('Enter the function\'s correct output: ', expected=float)
        print("Working...calculation may take some time depending on size of interval and desired accuracy")
        if settings is not None:
            # Correct func choice for 0 indexing
            result = iteration.shock_iterate(self.__gamma, func - 1, correct_output, settings[0], settings[1],
                                             settings[2], settings[3])
        else:
            result = iteration.shock_iterate(self.__gamma, func - 1, correct_output)
        if result is not None:
            print('M = %f' % result)
        else:
            print("No solution found within the specific interval")

    def expansion_iter_wrapper(self):
        """ Command line wrapper over iteration.expansion_iterate() """
        func = 1  # only expansion function is iteration.nu()
        settings = self.mbound_err_prompt()
        correct_output = prompt('Enter the correct value of nu: ', expected=float)
        print("Working...calculation may take some time depending on size of interval and desired accuracy")
        if settings is not None:
            # Correct func choice for 0 indexing
            result = iteration.expansion_iterate(self.__gamma, func - 1, correct_output, settings[0], settings[1],
                                                 settings[2], settings[3])
        else:
            result = iteration.expansion_iterate(self.__gamma, func - 1, correct_output)
        if result is not None:
            print('M = %f' % result)
        else:
            print("No solution found within the specific interval")

    @staticmethod
    def mbound_err_prompt():
        """
        Prompts user for mach number bound, and accuracy
        """
        need_args = confirm('Change iteration mach number range and accuracy (1<M<7, accuracy=1e-4) [y/N]? ')
        if need_args:
            start_M = prompt('Enter lower bound for mach #: ', expected=float)
            end_M = prompt('Enter upper bound for mach #: ', expected=float)
            accuracy = prompt('Enter accuracy (ex. 1e-4): ', expected=float)
            # Final element: step
            return [start_M, end_M, accuracy, accuracy * .1]
        else:
            pass


if __name__ == '__main__':
    CompressibleCLI()
