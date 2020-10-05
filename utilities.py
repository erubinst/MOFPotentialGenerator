# Common parameters
import numpy as np
from scipy.misc import derivative


# parameters in the electrostatic potential for H2
# units of angstroms squared times q
Q = 0.12
# k is used when using mks units - is in units of 1/(F/m)
k = 1/(4*np.pi*8.8*(10**-12))
# print(k)
# charge of electron in coulombs
qe = 1.6*10**-19
# alpha is in unites of Angstroms cubed
alpha = 0.675

# coordinates - Zn:0-2, O1:3, O2:4-9, C1:10-12
# gives distance - taken from rosnowCode
def magnitude(r,rspace):
    # 3D distance formula
    return np.sqrt((r[0]-rspace[0])**2+(r[1]-rspace[1])**2+(r[2]-rspace[2])**2)

def beta(C,T):
    return 1/(C*T)
# kB units in J*K^-1
kB = (1.38*10**-23)
# in units of J*s
hbar = (1.054*10**-34)

def deBroglieCoeff(isotope,T):
    return ((6*isotope.mass)/(beta(kB,T)*hbar**2))/(10**20)

def secondDeriv(f,x,dx):
    def firstDeriv(xx):
        return derivative(f,xx,dx)
    return derivative(firstDeriv, x, dx)


