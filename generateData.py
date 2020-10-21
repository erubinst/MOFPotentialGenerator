import atom as at
import realisticPotentialModel as pm
import copy as cp
import numpy as np
import time

def generatePotentialData(potentialType, xArray, yArray, zArray, atoms, source,quantump, T):
    start = time.perf_counter_ns()
    minimumy = []
    minimumx = []
    minimum = float("inf")
    potentials = []
    gaussian = []
    # loads in isotope object
    isotope = cp.copy(source)
    for val in zArray:
        # sets point of isotope to look at all points along z axis
        isotope.setPoint(0, 0, val)
        if quantump:
            u_sum = pm.UFH(potentialType, xArray, yArray, zArray, atoms, isotope, T)[0]
            gaussianPoint = pm.UFH(potentialType, xArray, yArray, zArray, atoms, isotope, T)[1]
            gaussian.append(gaussianPoint)
        else:
            # classsical version
            u_sum = pm.U(potentialType, isotope, atoms)
        if u_sum < minimum:
            minimum = u_sum
            zval = val
        potentials.append(u_sum)
    minimumy.append(minimum)
    minimumx.append(zval)
    end = time.perf_counter_ns()
    elapsedTime = (end-start)*(1/(10**9))
    print("Elapsed Time is", elapsedTime)
    return (potentials, gaussian)

#does not work in cases where the grid is not equal for all sides (|x|=|y|=|z|)
def saveFile(fileName, potentialType, xArray, yArray, zArray, atoms, source,quantump, T):
    potential = generatePotentialData(potentialType, xArray, yArray, zArray, atoms, source,quantump, T)[0]
    np.savetxt(fileName,np.transpose([xArray, yArray,zArray,potential]), delimiter = ", ")


saveFile("Potential22.csv", "Ulj", np.linspace(-3,3,50), np.linspace(-3,3,50), np.linspace(2,5,50), at.atoms, at.hydrogen, True, 22)
# saveFile("Potential77.csv","Ulj", np.linspace(-3,3,50), np.linspace(-3,3,50), np.linspace(2,5,50), at.atoms, at.hydrogen, True, 77)
