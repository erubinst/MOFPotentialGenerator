import matplotlib.pyplot as plt
import CuUtilities as ut

def plotPotential(potentialFile):
    dist = ut.readinPotential(potentialFile)[0]
    pot = ut.readinPotential(potentialFile)[1]
    minDist = ut.findMinimum(dist, pot)[0]
    minPot = ut.findMinimum(dist, pot)[1]
    plt.plot(dist, pot)
    plt.xlabel("Distance (angstroms)")
    plt.ylabel("Potential (K)")
    plt.text(0.02, .985, "Dist at Min (Angstroms):" + str(minDist), transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.text(0.02, .93, "Potential Min (K):" + str(minPot), transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.show()


def plotInterpolatedPotential(distFile, potFile):
    dist = ut.readinFile(distFile, False)
    pot = ut.readinFile(potFile,False)
    minDist = ut.findMinimum(dist, pot)[0]
    minPot = ut.findMinimum(dist, pot)[1]
    plt.plot(dist, pot)
    plt.xlabel("Distance (angstroms)")
    plt.ylabel("Potential (K)")
    plt.text(0.02, .985, "Dist at Min (Angstroms):" + str(minDist), transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.text(0.02, .93, "Potential Min (K):" + str(minPot), transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.show()


plotInterpolatedPotential("distance2.txt", "wB97.txt")
