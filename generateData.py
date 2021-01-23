import atom as at
import realisticPotentialModel as pm
import numpy as np

# does not work in cases where the grid is not equal for all sides (|x|=|y|=|z|)
def saveU2DFile(fileName, potentialType, xArray, yArray, zArray, atoms, source,quantump, T):
    potential = pm.generate2DPotentialData(potentialType, xArray, yArray, zArray, atoms, source,quantump, T)[0]
    np.savetxt(fileName, potential, delimiter=", ")

saveUFile("PotentialTotalH22.txt", "U", np.linspace(-3,3,50), np.linspace(-3,3,50), np.linspace(2,5,50), at.atoms, at.hydrogen, True, 22)
# saveUFile("PotentialTotal77.txt","U", np.linspace(-3,3,50), np.linspace(-3,3,50), np.linspace(2,5,50), at.atoms, at.deuterium, True, 77)
# saveUFile("PotentialClassicalFor22.txt", "U", np.linspace(-3,3,50), np.linspace(-3,3,50), np.linspace(3,7,50), at.atoms, at.deuterium, False,77)

# np.savetxt("22ZArray.txt",np.linspace(3,7,50), delimiter = ", ")

def saveU3DFile(fileName,potentialType, xArray, yArray, zArray, zval, atoms, source, quantump, T):
    potential = pm.generate3DPotentialData(potentialType, xArray, yArray, zArray, zval, atoms, source, quantump, T)[0]
    np.savetxt(fileName, potential,
    
