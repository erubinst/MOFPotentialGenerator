import atom as at
import realisticPotentialModel as pm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def UvZPlot(potentialTypes, atoms, T):
    z = np.linspace(2, 5)
    x = np.linspace(-3, 3)
    y = np.linspace(-3, 3)
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
        minText += minval[2] + " Minimum: " + str(round(minval[0],3)) + ", " + str(round(minval[1]))
        if index != lastMinIndex:
            minText += "\n"
    plt.text(0.02, .985, minText, transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.legend()
    plt.show()


def UvXPlot(potentialTypes, atoms, T,zmin):
    z = np.linspace(3, 5)
    x = np.linspace(-3, 3)
    y = np.linspace(-3, 3)
    minimums = []

    # Default size is about 6.5 x 4.75 (w x h)
    plt.figure(figsize=(8.125, 5.937))
    for pType in potentialTypes:

        def UxClassical(x):
            pots, xval_at_min, minPotential, elapsedTime = pm.generate2DXPotentialData(pType, x, y, z, atoms, at.hydrogen, False, T, zmin)
            minimums.append((xval_at_min, minPotential, "Classical " + pType))
            return pots

        def UxHFH(x):
            pots, xval_at_min, minPotential, elapsedTime = pm.generate2DXPotentialData(pType, x, y, z, atoms, at.hydrogen, True, T, zmin)
            minimums.append((xval_at_min, minPotential, "H2 " + pType))
            return pots

        def UxDFH(x):
            pots, xval_at_min, minPotential, elapsedTime = pm.generate2DXPotentialData(pType, x, y, z, atoms, at.deuterium, True, T,zmin)
            minimums.append((xval_at_min, minPotential, "D2 " + pType))
            return pots

        plt.plot(x, UxClassical(x), color="black", label="Classical " + pType)
        plt.plot(x, UxHFH(x), color="red", label="H2 " + pType)
        plt.plot(x, UxDFH(x), color="green", label="D2 " + pType)
    plt.xlabel("x axis distance (angstroms)")
    plt.ylabel("Potential (K)")
    lastMinIndex = len(minimums) - 1
    minText = ""
    for index in range(lastMinIndex + 1):
        minval = minimums[index]
        minText += minval[2] + " Minimum: " + str(round(minval[0],3)) + ", " + str(round(minval[1]))
        if index != lastMinIndex:
            minText += "\n"
    plt.text(0.02, .985, minText, transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.legend()
    plt.show()


def UvYPlot(potentialTypes, atoms, T,zmin):
    z = np.linspace(3, 5)
    x = np.linspace(-3, 3)
    y = np.linspace(-3, 3)
    minimums = []

    # Default size is about 6.5 x 4.75 (w x h)
    plt.figure(figsize=(8.125, 5.937))
    for pType in potentialTypes:

        def UyClassical(y):
            pots, yval_at_min, minPotential, elapsedTime = pm.generate2DYPotentialData(pType, x, y, z, atoms, at.hydrogen, False, T, zmin)
            minimums.append((yval_at_min, minPotential, "Classical " + pType))
            return pots

        def UyHFH(y):
            pots, yval_at_min, minPotential, elapsedTime = pm.generate2DYPotentialData(pType, x, y, z, atoms, at.hydrogen, True, T, zmin)
            minimums.append((yval_at_min, minPotential, "H2 " + pType))
            return pots

        def UyDFH(y):
            pots, yval_at_min, minPotential, elapsedTime = pm.generate2DYPotentialData(pType, x, y, z, atoms, at.deuterium, True, T,zmin)
            minimums.append((yval_at_min, minPotential, "D2 " + pType))
            return pots

        plt.plot(y, UyClassical(y), color="black", label="Classical " + pType)
        plt.plot(y, UyHFH(z), color="red", label="H2 " + pType)
        plt.plot(y, UyDFH(z), color="green", label="D2 " + pType)
    plt.xlabel("y axis distance (angstroms)")
    plt.ylabel("Potential (K)")
    lastMinIndex = len(minimums) - 1
    minText = ""
    for index in range(lastMinIndex + 1):
        minval = minimums[index]
        minText += minval[2] + " Minimum: " + str(round(minval[0],3)) + ", " + str(round(minval[1]))
        if index != lastMinIndex:
            minText += "\n"
    plt.text(0.02, .985, minText, transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor='red', alpha=0.5))
    plt.legend()
    plt.show()


def UvZPlot3D(potentialType, source, atoms, zval, quantump, T):
    z = np.linspace(2.5, 6, 100)
    x = np.linspace(-3, 3, 100)
    y = np.linspace(-3, 3, 100)
    X, Y = np.meshgrid(x, y)
    potentials, xval_at_minPotential, yval_at_minPotential, minPotential, elapsedTime = pm.generate3DPotentialData(potentialType, x, y, z, zval, atoms, source, quantump, T)
    # Default size is about 6.5 x 4.75 (w x h)
    fig = plt.figure(figsize=(8.125, 5.937))
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, np.array(potentials), rstride=1, cstride=1, cmap=cm.gist_rainbow)
    ax.set_xlabel("x axis distance (angstroms)")
    ax.set_ylabel("y axis distance (angstroms)")
    ax.set_zlabel(potentialType + " (K)")
    ax.text2D(0, 1, "U min: " + str(round(minPotential,3)) + "\nPoint: (" + str(round(xval_at_minPotential, 3)) + "," + str(round(yval_at_minPotential, 3)) + "," + str(round(zval,3)) + ")", transform=ax.transAxes, bbox=dict(facecolor='red', alpha=0.5))
    plt.show()


atoms = at.atoms
hydrogen = at.hydrogen
deuterium = at.deuterium

#UvZPlot(["U"], atoms, 22)
#UvZPlot3D("U", at.hydrogen, atoms, 3.04, False, 77)#UvXPlot(["U"], atoms, 22, 3.4)
UvYPlot(["U"], atoms, 22, 3.4)
