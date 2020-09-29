# Electric field computation
def correctionFactors(HAtom, AAtom, R):
    X = (HAtom.x-AAtom.x)
    Y = (HAtom.y-AAtom.y)
    Z = (HAtom.z-AAtom.z)

    xCorrection = X/R
    yCorrection = Y/R
    zCorrection = Z/R

    return(xCorrection, yCorrection, zCorrection)

def Efield(isotope, atom):
    mag = atom.magnitude(isotope)
    correctionFactor = correctionFactors(isotope, atom, mag)
    chargeFactor = atom.charge/(mag**2)
    Ex = chargeFactor*correctionFactor[0]
    Ey = chargeFactor*correctionFactor[1]
    Ez = chargeFactor*correctionFactor[2]
    # return (Ex, Ey, Ez)
    return (Ex**2+Ey**2+Ez**2)

# Different ways to compute potential
def Upol(isotope, atom):
    # return -((1/uT.kB)*(10**10)*uT.k*((1.6*10**-19)**2)*(0.675/2)*eF.Efield(isotope, r, charges, limit1, limit2))
    return -(0.675/2)*Efield(isotope, atom)

def Ulj(isotope, atom):
    R = atom.magnitude(isotope)
    epsilon = atom.jointEpsilon(isotope)
    sigma = atom.jointSigma(isotope)
    return 4*epsilon*((sigma/R)**12-(sigma/R)**6)
