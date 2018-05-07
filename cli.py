# CLI for compressible flow calculations
# By Ariel Mordoch, distributed under GNU GPL
import iteration
import util
import math

class compressibleCLI:

    def __init__(self):
        self.start()
        self.__gamma = 1.4

    def printRow(self, table, M):
        if table == 1:
            # Isentropic print
            A_As = util.AoverAstar(M, self.__gamma)
            T0_T = util.T0_over_T(M, self.__gamma)
            p0_p = util.p0_over_p(M, self.__gamma)
            rho0_rho = util.rho0_over_rho(M, self.__gamma)
            print("A/A* = %f \nT0/T = %f \np0/p = %f\nrho0/rho = %f" % (A_As, T0_T, p0_p, rho0_rho) )            
        elif table == 2 or table == 3:
            # Normal shocks/oblique shocks
            if table == 3:
                theta = float(input("Enter turning angle (deg): "))*math.pi/180
                alpha = input("For a strong oblique shock, enter 0: ")
                if alpha == "0": # Using string here avoids having a try-catch
                    beta = util.beta(M, self.__gamma, theta, 0)
                else:
                    beta = util.beta(M, self.__gamma, theta)
                M1n = M*math.sin(beta)
                M2 = util.M2(M1n, self.__gamma) # M2n
                M2_real = M2 / math.sin(beta - theta) # M2
            else:
                M2 = util.M2(M, self.__gamma)
                
            T2_T1 = util.T2overT1(M, M2, self.__gamma)
            rho2_rho1 = util.rho2_over_rho1(M, M2, self.__gamma)
            p2_p1 = util.p2_over_p1(M, M2, self.__gamma)
            p02_p01 = util.p02_over_p01(M, M2, self.__gamma)
            p02_p1 = util.p02_over_p1(M, M2, self.__gamma)
            rho02_rho01 = util.rho02_over_rho01(M, M2, self.__gamma)
            A2star_A1star = util.A2star_over_A1star(M, M2, self.__gamma)
            v2_v1 = util.v2overv1(M, M2, self.__gamma)
            if table == 3:
                # Some things must be modified for oblique shocks
                print("M1 = %f\nM1n = %f\nM2n = %f\nM2 = %f" % (M, M1n, M2, M2_real) )
                print("T2/T1 = %f\nrho2 / rho1 = %f\np2/p1 = %f\np02/p01 = %f" % (T2_T1, rho2_rho1, p2_p1, p02_p01) )
                print("rho02 / rho01 = %f\nA2*/A1* = %f\nv2/v1 = %f" % (rho02_rho01, A2star_A1star, v2_v1) )
            else:
                print("M2 = %f\nT2/T1 = %f\nrho2 / rho1 = %f\np2/p1 = %f\np02_p01 = %f" % (M2, T2_T1, rho2_rho1, p2_p1, p02_p01) )
                print("p02/p1 = %f\nrho02 / rho01 = %f\nA2*/A1* = %f\nv2/v1 = %f" % (p02_p1, rho02_rho01, A2star_A1star, v2_v1) )
        elif table == 4:
            # Expansions
            nu = util.nu(M, self.__gamma)*180/math.pi
            mu = util.mu(M)*180/math.pi
            print("nu = %f deg\nmu = %f deg" % (nu, mu) )
            
    def knownM(self):
        # Select isentropic, normal shock, oblique shock, or expansion
        print("Select option:")
        print("1: Isentropic Flow Properties\n2: Normal Shock Properties")
        print("3: Oblique Shock Properties\n4: Prandtl-Meyer\n5: Back")
        # Inputs & input validation
        try:
            # Several layers of input checking here. First, the initial selection
            sel = int(input())
        except Exception as e:
            # Restart if there was a mistake during selection
            print("Error: ", e)
            print("Continuing...")
            self.knownM()
        # Input validation for mach number  
        try:
            if sel != 5:
                M = float(input("Enter Mach number: "))
                # Use this variable to avoid needing nested try-catch statements
                # in case an exception occurs later/during recursion
                returnToPrev = False
            else:
                returnToPrev = True
        except Exception as e:
                print("Error: ", e)
                print("Continuing...")
                self.knownM()
        else:
            if not(returnToPrev):
                # Additional error handling required for mathematical exceptions
                try:
                    self.printRow(sel, M)
                except Exception as e:
                    print("Error: ", e)
                    print("Continuing...")

            
                
    def start(self):
        # Intro
        print("Compressible Flow Calculator 1.0")
        print("By Ariel Mordoch")
        try:
            self.__gamma = float(input("Enter a value for gamma (leave blank for 1.4): "))
        except ValueError:
            print("Assuming gamma = 1.4")
            self.__gamma = 1.4
        # Start Looping
        loop = True
        while loop:
            # Decide between iteration & calculations
            print("---\nSelect option: ")
            print("1: Known Mach number")
            print("2: Unknown Mach number")
            print("3: Exit")
            try:
                sel = int(input())
            except ValueError:
                print("Goodbye!")
                loop = False
            else:
                if sel == 1:
                    self.knownM()
                elif sel == 2:
                    self.unknownM()
                else:
                    print("Goodbye!")
                    loop = False
                
            
            
CLI = compressibleCLI() 
            
        
