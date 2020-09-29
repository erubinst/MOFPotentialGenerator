import numpy as np

class point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def setPoint(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def distance(self, pt):
        return np.sqrt((self.x-pt.x)**2+(self.y-pt.y)**2+(self.z-pt.z)**2)

def convertTupleToPoint(tuple):
    return point(tuple[0], tuple[1], tuple[2])
