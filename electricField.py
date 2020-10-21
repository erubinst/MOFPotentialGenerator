# Esme Rubinstein
# Finds electric field in each direction for a series of point charges

import utilities as uT
import correctionFactors as cF

def Efield(point, r, charges, limit1, limit2):
    Ex = 0
    Ey = 0
    Ez = 0
    for i in range(limit1, limit2+1):
        mag = uT.magnitude(point, r[i])
        correctionFactor = cF.correctionFactors(point, r[i], mag)
        chargeFactor = charges[i]/(mag**2)
        Ex += chargeFactor*correctionFactor[0]
        Ey += chargeFactor*correctionFactor[1]
        Ez += chargeFactor*correctionFactor[2]
    # return (Ex, Ey, Ez)
    return (Ex**2+Ey**2+Ez**2)

R = [3.39563283, 3.39563283, 3.39563283, 3.5, 3.230114895, 3.230114895, 3.230114895, 3.230114895, 3.230114895, 3.230114895, 3.255457303, 3.255457303, 3.255457303]
