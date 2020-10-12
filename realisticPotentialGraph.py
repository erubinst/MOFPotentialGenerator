import atom as at
import realisticPotentialModel as pm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import copy as cp
import time


def UvZPlot(potentialType, atoms, T):
    z = np.linspace(3, 5, 30)
    x = np.linspace(-3, 3, 30)
    y = np.linspace(-3, 3, 30)
    minimums = []

    # Source is for choosing hydrogen or deuterium and atoms is the list of atoms in MOF
    def UzInternal(z, quantump, source):
        minimum = float("inf")
        potentials = []
        # loads in isotope object
        isotope = cp.copy(source)
        for val in z:
            # sets point of isotope to look at all points along z axis
            isotope.setPoint(0, 0, val)
            if quantump:
                u_sum = pm.UFH(potentialType, x, y, z, atoms, isotope, T)
            else:
                # classsical version
                u_sum = pm.U(potentialType, isotope, atoms)
            if u_sum < minimum:
                minimum = u_sum
                zval = val
            potentials.append(u_sum)
        print("Minimum Point: (", zval, ",", int(minimum), ")")
        minimums.append(minimum)
        return potentials

    def UzClassical(z):
        return UzInternal(z, False, at.hydrogen)

    def UzHFH(z):
        return UzInternal(z, True, at.hydrogen)

    def UzDFH(z):
        return UzInternal(z, True, at.deuterium)

    start = time.perf_counter()
    # Default size is about 6.5 x 4.75 (w x h)
    plt.figure(figsize=(8.125,5.937))
    plt.plot(z, UzClassical(z), label="Classical")
    plt.plot(z, UzHFH(z), label="H2")
    plt.plot(z, UzDFH(z), label="D2")
    end = time.perf_counter()
    print("Elapsed Time is", end - start)
    plt.xlabel("z axis distance (angstroms)")
    plt.ylabel(potentialType + " (K)")
    plt.text(0.02, .985, "Classical Minimum: " + str(round(minimums[0])) + "\nH2 Minimum: " + str(round(minimums[1])) + "\nD2 Minimum: " + str(round(minimums[2])), transform=plt.gcf().transFigure, va='top',bbox=dict(facecolor='red', alpha=0.5))
    plt.legend()
    plt.show()


def UvZPlot3D(potentialType, source, atoms, zval, T):
    z = np.linspace(2.5, 6, 20)
    x = np.linspace(-3, 3, 20)
    y = np.linspace(-3, 3, 20)

    def Uz3DInternal(x, y, z, quantump):
        potentials = []
        minimum = float("inf")
        yMinVal = 0
        xMinVal = 0
        isotope = cp.copy(source)
        for yval in y:
            row_potentials = []
            for xval in x:
                isotope.setPoint(xval, yval, zval)
                if quantump:
                    u_sum = pm.UFH(potentialType, x, y, z, atoms, isotope, T)
                else:
                    # classical version
                    u_sum = pm.U(potentialType, isotope, atoms)
                u_sum = min(u_sum, 4000)
                if u_sum < minimum:
                    minimum = u_sum
                    yMinVal = yval
                    xMinVal = xval
                row_potentials.append(u_sum)
            potentials.append(row_potentials)
        return np.array((np.array(potentials),minimum,xMinVal,yMinVal), dtype=object)

    X, Y = np.meshgrid(x, y)
    start = time.perf_counter()
    Z = Uz3DInternal(x, y, z, True)
    end = time.perf_counter()
    print("Elapsed Time is", end - start)
    # Default size is about 6.5 x 4.75 (w x h)
    fig = plt.figure(figsize=(8.125,5.937))
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z[0], rstride=1, cstride=1, cmap=cm.gist_rainbow)
    ax.set_xlabel("x axis distance (angstroms)")
    ax.set_ylabel("y axis distance (angstroms)")
    ax.set_zlabel(potentialType + " (K)")
    ax.text2D(0, 1, "U min: " + str(int(Z[1])) + "\nPoint: (" + str(round(Z[2], 3)) + "," + str(round(Z[3], 3)) + "," + str(round(zval,2)) + ")", transform=ax.transAxes, bbox=dict(facecolor='red', alpha=0.5))
    plt.show()


atoms = at.atoms
hydrogen = at.hydrogen
deuterium = at.deuterium

UvZPlot("Ulj", atoms, 77)
# UvZPlot3D("Ulj", at.hydrogen, atoms, 3.3, 77)
