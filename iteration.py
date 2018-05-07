# Compressible flow calculator - iteration functions
# By Ariel Mordoch - distributed under GNU GPL3
from util import *
import math
from numpy import linspace


def iterate(correctOutput, x, y, accuracy):
    """
    -- Description --
    iterate() takes in a given function output [y] = [f(x)],
    an array of input values [x], and an array of output values [f(x)],
    and iterates through the array to find
    where the desired output occurs.
    -- Parameters --
    correctOutput - desired function output
    x - array of input values
    y - array of output values
    accuracy - desired percent error
    """
    accuracy = abs(accuracy)
    i = 0;
    foundX = None;
    for output in y:
        pctErr = abs( (correctOutput - output)/correctOutput )
        if pctErr <= accuracy:
            foundX = x[i]
            break
        else:
            i+=1
    return foundX

# Iteration over isentropic functions #

def isentropicIterate(gamma, correctOutput, startM, endM, accuracy=1e-4, step=1e-5):
    # Iterates over isentropic flow functions using input step size and accuracy
    # Create array of inputs
    machNums = []
    for i in linspace(startM, endM, (endM - startM)/step + 1):
        machNums.append(i)
    # Select function to iterate over
    print("Select function to iterate over:")
    print("1: A / A*")
    print("2: T0 / T")
    print("3: p0 / p")
    print("4: rho0 / rho")
    func = int( input("?: ") )
    print("Working...calculation may take some time depending on size of interval and desired accuracy")
    outputs = []
    # Generate array of outputs
    if func == 1:
        for M in machNums:
            outputs.append( AoverAstar(M, gamma) )
    elif func == 2:
        for M in machNums:
            outputs.append( T0_over_T(M, gamma) )
    elif func == 3:
        for M in machNums:
            
            outputs.append( p0_over_p(M, gamma) )
    elif func == 4:
        for M in machNums:
            outputs.append( rho0_over_rho(M, gamma) )
    # Finally, iterate.
    foundM = iterate(correctOutput, machNums, outputs, accuracy)
    if foundM == None:
        print("No solution found within the specified interval")
    else:
        return round(foundM, len(str(accuracy)) )

# Iteration over normal shock functions #

def shockIterate(gamma, correctOutput, startM, endM, accuracy=1e-4, step=1e-5):
    # Iterates over normal shock functions using input step size and accuracy
    # Create array of inputs
    machNums = []
    for i in linspace(startM, endM, (endM - startM)/step + 1):
        machNums.append(i)
    # Select function to iterate over
    print("Select function to iterate over:")
    print("1: p02 / p01")
    print("2: T2 / T1")
    print("3: p2 / p1")
    func = int( input("?: ") )
    print("Working...calculation may take some time depending on size of interval and desired accuracy")
    outputs = []
    # Generate array of outputs
    if func == 1:
        for M1 in machNums:
            outputs.append( p02_over_p01(M1 , M2(M1, gamma), gamma) ) 
    elif func == 2:
        for M1 in machNums:
            outputs.append( T2overT1(M1, M2(M1, gamma), gamma) )
    elif func == 3:
        for M1 in machNums:
            outputs.append( p2_over_p1(M1, M2(M1, gamma), gamma) )
    # Finally, iterate.
    foundM = iterate(correctOutput, machNums, outputs, accuracy)    
    if foundM == None:
        print("No solution found within the specified interval")
    else:
        return round(foundM, len(str(accuracy)) )

# Iteration over expansion functions #
def expansionIterate(gamma, correctOutput, startM, endM, accuracy=1e-4, step=1e-5):
    # Iterates over nu(M, gamma) to find the mach number
    machNums = []
    for i in linspace(startM, endM, (endM-startM)/step+1):
        machNums.append(i)
    # Generate array of nu values
    nuVals = []
    for M in machNums:
        nuVals.append(nu(M, gamma))
    foundM = iterate(correctOutput, machNums, nuVals, accuracy)
    if foundM == None:
         print("No solution found within the specified interval")
    else:
        return round(foundM, len(str(accuracy)) )

        
        
