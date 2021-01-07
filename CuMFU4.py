import matplotlib.pyplot as plt

def readinPotential(file):
    file1 = open(file, "r")
    distance = []
    potential = []
    numlist = file1.read().split()
    for i in range(len(numlist)):
        numlist[i] = float(numlist[i])
        if i%2 == 0:
            distance.append(numlist[i])
        else:
            numlist[i] = numlist[i]/.008314
            potential.append(numlist[i])
    file1.close()
    return distance,potential

def readinFile(file, potp):
    file1 = open(file, "r")
    newlst = []
    numList = file1.read().split()
    for i in range(len(numList)):
        numList[i] = float(numList[i])
        if potp == True:
            numList[i] = numList[i]/.008314
        newlst.append(numList[i])
    file1.close()
    return newlst

#print(readinFile("distance01.txt", False))

def findMinimum(dist, pot):
    minimumPot = 100000
    minimumDist = 0
    for i in range(len(pot)):
        if pot[i] < minimumPot:
            minimumPot = pot[i]
            minimumDist = dist[i]
    return minimumDist, minimumPot


def plotPotential(potentialFile):
    dist = readinPotential(potentialFile)[0]
    pot = readinPotential(potentialFile)[1]
    minDist = findMinimum(dist, pot)[0]
    minPot = findMinimum(dist, pot)[1]
    plt.plot(dist, pot)
    plt.xlabel("Distance (angstroms)")
    plt.ylabel("Potential (K)")
    plt.text(0.02, .985, "Dist at Min (Angstroms):" + str(minDist), transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.text(0.02, .93, "Potential Min (K):" + str(minPot), transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.show()

def plotInterpolatedPotential(distFile, potFile):
    dist = readinFile(distFile, False)
    pot = readinFile(potFile,True)
    minDist = findMinimum(dist, pot)[0]
    minPot = findMinimum(dist, pot)[1]
    plt.plot(dist, pot)
    plt.xlabel("Distance (angstroms)")
    plt.ylabel("Potential (K)")
    plt.text(0.02, .985, "Dist at Min (Angstroms):" + str(minDist), transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.text(0.02, .93, "Potential Min (K):" + str(minPot), transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.show()


probGS = readinFile("probGS.txt", False)
probFE = readinFile("probFE.txt", False)

areaGS = sum(probGS)
areaFE = sum(probFE)

print(areaGS)
