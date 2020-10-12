# Common constants and functions for computing potential energy.

# Esme Rubinstein

from scipy.misc import derivative
# import scipy as sym
import numpy as np

# distance points??
# r = np.linspace(-5*10**-10,5*10**-10)
r = np.linspace(2.8,7,100)

# parameters
# sigma is in Angstroms
sigma = 2.8445
# epsilon is in kJ/mol - i know this is depth of the well, found 35 in MFU breakthrough paper but this could be wrong
epsilon = 25.85
# prefactor parts in mks and then converted at the end
# reduced mass - can alter to look at D2 or H2
muH = 2*(1.7*10**-27)
muD = 2*muH
# T = 600
def beta(C,T):
    return 1/(C*T)


# gas constant R in units kJ/mol*K
R = 8.314*0.001
# kB units in J*K^-1
kB = (1.38*10**-23)
# in units of J*s
hbar = (1.054*10**-34)
# conversion factor of 10**20 to go from m**2 to A**2
def prefac1(mu, T):
    return ((beta(kB,T)*hbar**2)/(24*mu))*(10**20)
def prefac2(mu, T):
    return (((beta(kB,T)**2)*(hbar**4))/(1152*(mu**2)))*(10**40)

# LJ Potential
def ULJ(r,epsilon,sigma):
    return 4*epsilon*((sigma/r)**12-(sigma/r)**6)

# derivative check - manually doing the derivative because I think my python derivs are wrong
def firstDerivCheck(r):
    return 4*epsilon*(6*(sigma**6))/(r**7)-(12*(sigma**12)/(r**13))
# manual second derivative just to get results
def secondDerivCheck(r):
    return -(168*epsilon*(sigma**6)*(r**6)-624*epsilon*(sigma**12))/(r**14)

def thirdDerivCheck(r):
    return (1344*epsilon*(sigma**6)*(r**6)-8736*epsilon*(sigma**12))/(r**15)

def fourthDerivCheck(r):
    return -(12096*epsilon*(sigma**6)*(r**6)-131040*epsilon*(sigma**12))/(r**16)

# derivative function - doesn't work - I'm guessing my function has some kind of aspect that doesn't work with python's built in library
def d(x):
    return derivative(ULJ, x)

#Potential Energy using Correction factors, set mu for molecule type
def Correction(mu,r,T):
    return ULJ(r,epsilon,sigma) + prefac1(mu,T)*((secondDerivCheck(r)+2*firstDerivCheck(r)/r)+prefac2(mu,T)*(15*firstDerivCheck(r)/(r**3)+(4*thirdDerivCheck(r)/r)+fourthDerivCheck(r)))

# print(prefac1(muH,77)*firstDerivCheck(6)/6)
