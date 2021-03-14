import copy as cp
import numpy as np
import utilities as uT
import time
import atom as at

# Electric field computation
def correctionFactors(HAtom, AAtom, R):
    X = (HAtom.x-AAtom.x)
    #print("x", X)
    #print("R", R)
    Y = (HAtom.y-AAtom.y)
    Z = (HAtom.z-AAtom.z)

    xCorrection = X/R
    yCorrection = Y/R
    zCorrection = Z/R

    return(xCorrection, yCorrection, zCorrection)


def Efield(isotope, atoms):
    Ex = 0
    Ey = 0
    Ez = 0
    for atom in atoms:
        R = atom.magnitude(isotope)
        correctionFactor = correctionFactors(isotope, atom, R)
        chargeFactor = atom.charge/(R**2)
        Ex += chargeFactor*correctionFactor[0]
        Ey += chargeFactor*correctionFactor[1]
        Ez += chargeFactor*correctionFactor[2]
    # return (Ex, Ey, Ez)
    # print("efield",Ex**2+Ey**2+Ez**2)
    # print("R", R)
    return (Ex**2+Ey**2+Ez**2)


# Different ways to compute potential
def Upol(isotope, atoms):
    # return -((1/uT.kB)*(10**10)*uT.k*((1.6*10**-19)**2)*(0.675/2)*eF.Efield(isotope, r, charges, limit1, limit2))
    # print(-((1/uT.kB)*(10**10)*uT.k*((uT.qe)**2)*(0.675/2)*Efield(isotope, atoms)))
    # print("upol", -((1/uT.kB)*(10**10)*uT.k*((uT.qe)**2)*(0.675/2)*Efield(isotope, atoms)))
    return (-((1/uT.kB)*(10**10)*uT.k*((uT.qe)**2)*(0.675/2)*Efield(isotope, atoms)))


def Ulj(isotope, atoms):
    Ulj = 0
    for atom in atoms:
        R = atom.magnitude(isotope)
        epsilon = atom.jointEpsilon(isotope)
        sigma = atom.jointSigma(isotope)
        Ulj += 4*epsilon*((sigma/R)**12-(sigma/R)**6)
    # print("ulj", Ulj)
    return Ulj


def U(potentialType, isotope, atoms):
    if potentialType == "Upol":
        return Upol(isotope, atoms)
    elif potentialType == "Ulj":
        return Ulj(isotope, atoms)
    elif potentialType == "U":
        #print("fullU", Ulj(isotope, atoms) + Upol(isotope, atoms))
        return Ulj(isotope, atoms) + Upol(isotope, atoms)


# computes the average of the classical potential weighted by a gaussian
def UFH(potentialType, xArray, yArray, zArray, atoms, source, T):
    u_sum = 0
    normalizationFactor = 0
    isotope = cp.copy(source)
    for gridZ in zArray:
        for gridY in yArray:
            for gridX in xArray:
                isotope.setPoint(gridX, gridY, gridZ)
                u_classical = U(potentialType, isotope, atoms)
                decayFactor = np.exp(-uT.deBroglieCoeff(isotope,T)*isotope.magnitude(source))
                u_sum += u_classical*decayFactor
                normalizationFactor += decayFactor
    return (u_sum/normalizationFactor)

# def UFH(potentialType, xArray, yArray, zArray, atoms, source, T):
#     u_sum = 0
#     normalizationFactor = 0
#     isotope = cp.copy(source)
#     for x, y, z in it.product(xArray, yArray, zArray):
#         isotope.setPoint(x, y, z)
#         u_classical = U(potentialType, isotope, atoms)
#         decayFactor = np.exp(-uT.deBroglieCoeff(isotope,T)*isotope.magnitude(source))
#         u_sum += u_classical*decayFactor
#         normalizationFactor += decayFactor
#     return u_sum/normalizationFactor

def generate2DPotentialData(potentialType, xArray, yArray, zArray, atoms, source, quantump, T):
    start = time.perf_counter()
    minPotential = float("inf")
    zval_at_minPotential = 0
    potentials = []
    # loads in isotope object
    isotope = cp.copy(source)
    for val in zArray:
        print(val)
        # sets point of isotope to look at all points along z axis
        isotope.setPoint(0, 0, val)
        if quantump:
            u_sum =UFH(potentialType, xArray, yArray, zArray, atoms, isotope, T)
        else:
            # classsical version
            u_sum = U(potentialType, isotope, atoms)
        if u_sum < minPotential:
            minPotential = u_sum
            zval_at_minPotential = val
        potentials.append(u_sum)
    end = time.perf_counter()
    elapsedTime = (end-start)*(1/(10**3))
    return (potentials, zval_at_minPotential, minPotential, elapsedTime)

def generate2DXPotentialData(potentialType, xArray, yArray, zArray, atoms, source, quantump, T,zval_at_min):
    start = time.perf_counter()
    minPotential = float("inf")
    xval_at_minPotential = 0
    potentials = []
    # loads in isotope object
    isotope = cp.copy(source)
    for val in xArray:
        # sets point of isotope to look at all points along z axis
        isotope.setPoint(val, 0, zval_at_min)
        if quantump:
            u_sum =UFH(potentialType, xArray, yArray, zArray, atoms, isotope, T)
        else:
            # classsical version
            u_sum = U(potentialType, isotope, atoms)
        if u_sum < minPotential:
            minPotential = u_sum
            xval_at_minPotential = val
        potentials.append(u_sum)
    end = time.perf_counter()
    elapsedTime = (end-start)*(1/(10**3))
    return (potentials, xval_at_minPotential, minPotential, elapsedTime)

def generate2DYPotentialData(potentialType, xArray, yArray, zArray, atoms, source, quantump, T,zval_at_min):
    start = time.perf_counter()
    minPotential = float("inf")
    yval_at_minPotential = 0
    potentials = []
    # loads in isotope object
    isotope = cp.copy(source)
    for val in yArray:
        # sets point of isotope to look at all points along z axis
        isotope.setPoint(0, val, zval_at_min)
        if quantump:
            u_sum =UFH(potentialType, xArray, yArray, zArray, atoms, isotope, T)
        else:
            # classsical version
            u_sum = U(potentialType, isotope, atoms)
        if u_sum < minPotential:
            minPotential = u_sum
            yval_at_minPotential = val
        potentials.append(u_sum)
    end = time.perf_counter()
    elapsedTime = (end-start)*(1/(10**3))
    return (potentials, yval_at_minPotential, minPotential, elapsedTime)

def generate3DPotentialData(potentialType, xArray, yArray, zArray, zval, atoms, source, quantump, T):
    start = time.perf_counter()
    potentials = []
    minPotential = float("inf")
    yval_at_minPotential = 0
    xval_at_minPotential = 0
    isotope = cp.copy(source)
    for yval in yArray:
        row_potentials = []
        for xval in xArray:
            isotope.setPoint(xval, yval, zval)
            if quantump:
                u_sum = UFH(potentialType, xArray, yArray, zArray, atoms, isotope, T)
            else:
                # classical version
                u_sum = U(potentialType, isotope, atoms)
            u_sum = min(u_sum, 4000)
            if u_sum < minPotential:
                minPotential = u_sum
                yval_at_minPotential = yval
                xval_at_minPotential = xval
            row_potentials.append(u_sum)
        potentials.append(row_potentials)
    end = time.perf_counter()
    elapsedTime = (end-start)*(1/(10**3))
    return (potentials, xval_at_minPotential, yval_at_minPotential, minPotential, elapsedTime)


def UzHydrogen(z):
    isotope = at.deuterium
    isotope.setPoint(0,0,z)
    return UFH("U", np.linspace(-3,3), np.linspace(-3,3), np.linspace(3,4), at.atoms, isotope, 22)

#print(generate2DPotentialData("U", np.linspace(-3,3,100), np.linspace(-3,3,100), np.linspace(2.4,6,150), at.atoms, at.hydrogen, False, 22))

#print(generate2DXPotentialData("U", np.linspace(-3,3), np.linspace(-3,3), np.linspace(2.4,6), at.atoms, at.hydrogen, False, 22, 3.04))

#print(np.linspace(2.4,6,150))
