import atom as at
import realisticPotentialModel as pm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import copy as cp
import time


def UvZPlot(potentialTypes, atoms, T):
    z = np.linspace(2, 7, 10)
    x = np.linspace(-3, 3, 10)
    y = np.linspace(-3, 3, 10)
    minimums = []

    # Default size is about 6.5 x 4.75 (w x h)
    plt.figure(figsize=(8.125, 5.937))
    for pType in potentialTypes:

        def UzClassical(z):
            pots, zval_at_min, minPotential, elapsedTime = pm.generate2DPotentialData(pType, x, y, z, atoms, at.hydrogen, False, T)
            minimums.append((zval_at_min, minPotential, "Classical " + pType))
            return pots

        def UzHFH(z):
            pots, zval_at_min, minPotential, elapsedTime = pm.generate2DPotentialData(pType, x, y, z, atoms, at.hydrogen, True, T)
            minimums.append((zval_at_min, minPotential, "H2 " + pType))
            return pots

        def UzDFH(z):
            pots, zval_at_min, minPotential, elapsedTime = pm.generate2DPotentialData(pType, x, y, z, atoms, at.deuterium, True, T)
            minimums.append((zval_at_min, minPotential, "D2 " + pType))
            return pots

        plt.plot(z, UzClassical(z), color="black", label="Classical " + pType)
        plt.plot(z, UzHFH(z), color="red", label="H2 " + pType)
        plt.plot(z, UzDFH(z), color="green", label="D2 " + pType)
    plt.xlabel("z axis distance (angstroms)")
    plt.ylabel("Potential (K)")
    lastMinIndex = len(minimums) - 1
    minText = ""
    for index in range(lastMinIndex + 1):
        minval = minimums[index]
        minText += minval[2] + " Minimum: " + str(round(minval[0])) + ", " + str(round(minval[1]))
        if index != lastMinIndex:
            minText += "\n"
    plt.text(0.02, .985, minText, transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.legend()
    plt.show()


def UvZPlot3D(potentialType, source, atoms, zval, quantump, T):
    z = np.linspace(2.5, 6, 5)
    x = np.linspace(-3, 3, 5)
    y = np.linspace(-3, 3, 5)
    X, Y = np.meshgrid(x, y)
    potentials, xval_at_minPotential, yval_at_minPotential, minPotential, elapsedTime = pm.generate3DPotentialData(potentialType, x, y, z, zval, atoms, source, quantump, T)
    # Default size is about 6.5 x 4.75 (w x h)
    fig = plt.figure(figsize=(8.125, 5.937))
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, np.array(potentials), rstride=1, cstride=1, cmap=cm.gist_rainbow)
    ax.set_xlabel("x axis distance (angstroms)")
    ax.set_ylabel("y axis distance (angstroms)")
    ax.set_zlabel(potentialType + " (K)")
    ax.text2D(0, 1, "U min: " + str(int(xval_at_minPotential)) + "\nPoint: (" + str(round(yval_at_minPotential, 3)) + "," + str(round(zval, 2)) + "," + str(round(minPotential,3)) + ")", transform=ax.transAxes, bbox=dict(facecolor='red', alpha=0.5))
    plt.show()


atoms = at.atoms
hydrogen = at.hydrogen
deuterium = at.deuterium

UvZPlot(["U","Ulj"], atoms, 77)
# UvZPlot3D("Ulj", at.hydrogen, atoms, 3.3, True, 77)
