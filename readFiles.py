def readinFile(file):
    file1 = open(file, "r")
    newlst = []
    numList = file1.read().split()
    for i in range(len(numList)):
        numList[i] = float(numList[i])
        newlst.append(numList[i])
    file1.close()
    return newlst

