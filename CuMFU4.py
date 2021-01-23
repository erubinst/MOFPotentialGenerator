import SE as se
import CuUtilities as ct
import atom as at

def generateWavefunction(potentialFile, distanceFile, isotope, GSE, FEE):
    potential = ct.readinFile(potentialFile, False)
    distance = ct.readinFile(distanceFile, False)
    psiGS = se.SE(potential, 0, 0.001, distance, isotope.mass, GSE)
    psiFE = se.SE(potential, 0, 0.001, distance, isotope.mass, FEE)
    probGS = ct.squaredlist(psiGS)
    probFE = ct.squaredlist(psiFE)
    areaGS = sum(probGS[0:184])
    areaFE = sum(probFE[0:184])
    #print(areaGS)
    #print(areaFE)
    ct.normLst(psiGS, areaGS)
    ct.normLst(psiFE, areaFE)
    ct.normLst(probGS, areaGS)
    ct.normLst(probFE, areaFE)
    #print(sum(probGS))
    #print(sum(probFE))
    return (psiGS, probGS, psiFE, probFE)

def main():
    #psiGSH,probGSH, psiFEH, probFEH = generateWavefunction("potential01.txt", "distance01.txt", at.deuterium, -2312.186611170497, -1100.063)
    psiGSH,probGSH, psiFEH, probFEH = generateWavefunction("B3lyp.txt", "distance2.txt", at.deuterium ,-4592.6687261555,-3495.8614957392012457)
    ct.makeFile("B3psiGSD.txt", psiGSH)
    ct.makeFile("B3probGSD.txt", probGSH)
    ct.makeFile("B3psiFED.txt", psiFEH)
    ct.makeFile("B3probFED.txt", probFEH)

#main()

def bumpPercentage(potentialFile, distanceFile, isotope, GSE, FEE, startPoint, endPoint):
    probFE = generateWavefunction(potentialFile, distanceFile, isotope, GSE, FEE)[3]
    bumpPercent = ct.sumLstSection(probFE, startPoint, endPoint)
    return bumpPercent

#print(bumpPercentage("potential01.txt", "distance01.txt", at.hydrogen, -2101.47875, -907.6, 0, 30))
#print(bumpPercentage("potential01.txt", "distance01.txt", at.hydrogen, -2101.47875, -907.6, 31, 74))
#print(bumpPercentage("potential01.txt", "distance01.txt", at.hydrogen, -2101.47875, -907.6, 75, 123))
#print(bumpPercentage("potential01.txt", "distance01.txt", at.hydrogen, -2101.47875, -907.6, 124, 213))


#print(bumpPercentage("potential01.txt", "distance01.txt", at.deuterium, -2312.18661117049, -1100.063, 0, 25))
#print(bumpPercentage("potential01.txt", "distance01.txt", at.deuterium, -2312.18661117049, -1100.063, 26, 48))
#print(bumpPercentage("potential01.txt", "distance01.txt", at.deuterium, -2312.18661117049, -1100.063, 49, 114))
#print(bumpPercentage("potential01.txt", "distance01.txt", at.deuterium, -2312.18661117049, -1100.063, 115, 213))

#print(bumpPercentage("B3lyp.txt", "distance2.txt", at.hydrogen, -4350.08496862226, -2903.24140354, 0, 43))
#print(bumpPercentage("B3lyp.txt", "distance2.txt", at.hydrogen, -4350.08496862226, -2903.24140354, 43, 228))

print(bumpPercentage("B3lyp.txt", "distance2.txt", at.deuterium, -4592.6687261555, -3495.8614957392, 0, 40))

print(bumpPercentage("B3lyp.txt", "distance2.txt", at.deuterium, -4592.6687261555, -3495.8614957392,41, 183))
