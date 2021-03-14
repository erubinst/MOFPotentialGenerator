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

    def rotate(self, rotAxis, angle):
        angle *= np.pi/180
        if rotAxis == 'x':
            originaly = self.y
            self.y = self.y*np.cos(angle)+self.z*np.sin(angle)
            self.z = -originaly*np.sin(angle)+self.z*np.cos(angle)
        if rotAxis == 'y':
            originalx = self.x
            self.x = self.x*np.cos(angle)-self.z*np.sin(angle)
            self.z = originalx*np.sin(angle)+self.z*np.cos(angle)

    def translation(self, axis, amount):
        if axis == 'x':
            self.x += amount
        if axis == 'y':
            self.y += amount
        if axis == 'z':
            self.z += amount

    def flipDirection(self):
        self.x = -self.x
        self.y = -self.y
        self.z = -self.z


def convertTupleToPoint(tuple):
    return point(tuple[0], tuple[1], tuple[2])
