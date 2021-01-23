# utility functions for Cu Potential Program
import numpy as np

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

def findMinimum(dist, pot):
    minimumPot = 100000
    minimumDist = 0
    for i in range(len(pot)):
        if pot[i] < minimumPot:
            minimumPot = pot[i]
            minimumDist = dist[i]
    return minimumDist, minimumPot

def sumLstSection(lst, start, end):
    sum = 0
    for i in range(start, end+1):
        sum += lst[i]
    return sum

def squaredlist(lst):
    lstSquared = []
    for i in range(len(lst)):
        lstSquared.append(lst[i]**2)
    return lstSquared


def makeFile(fileName, Lst):
    np.savetxt(fileName, Lst, delimiter=", ")


def normLst(lst, area):
    for i in range(len(lst)):
        lst[i] = lst[i]/area


