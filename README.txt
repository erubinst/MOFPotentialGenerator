# MOFPotentialGenerator
This is the code that is a part of my Honors research at Oberlin College
CuMFU - returns normalized wavefunction psi in a text file for given potential file for the CuMFU, calculates area under psi^2 for a given range
CuUtilities - contains utility functions specific to CuMFU potential files, such as reading in potential text file into a list
SE - Utilizes wag the dog method to find energy levels for ground state and first excited state for any general MOF or MFU.  plugs a given potential into the schrodinger equation, potential must have equal distance spacing
atom - class definition of an atom in a MOF (MOF 5 is defined in the file below the atom class definition), sigma and epsilon values for MOF5 taken from experimental parameters in 2009 paper from UC Berkeley
electricField - old file used to calculate summed electric field for calculating potential effects from polarizability.  file not usable since we discovered that the FH approximation must be applied to potential before summing.
generateData - saves text files of generated potentials.  used mostly for MOF5
PlotCuPot - plots potential file for CuMFU given by UC Berkeley, calculates minimum for file
point - defines point class.  used for modelling a given MOF and is used in atom class
potentialEnergy - preliminary file when I was first starting out.  models a simple 1D lennard-jones potential and runs potential through FH equation for H2 and D2
readFiles - preliminary readFiles function, not used currently
realisticPotentialGraph - graphing functions for 1D and 3D potentials for a given MOF (used mostly for MOF5 where we have parameters). 
realisticPotentialModel - generates a model of the potential for a given MOF (a MOF is inputed as a series of atom class objects) using lennard jones and polarizability energy plugged into a version of the FH that averages over a Gaussian wavefunction
rosnowPotentialGenerator - code from previous student with other modeling techniques for MOF potential, contains coordinates used for MOF5
utilities - contains common constants and a function for calculating second derivative

