from decimal import ROUND_HALF_DOWN
from sklearn.feature_selection import SequentialFeatureSelector

import numpy

class Mesh:
    def __init__(self,width,height,c,nu,rho,dt):
        self.dx = 2 / (width - 1)
        self.dy = 2 / (height - 1)

        x = numpy.linspace(0, 2, width)
        y = numpy.linspace(0, 2, height)
        self.X, self.Y = numpy.meshgrid(x, y)

        self.u = numpy.zeros((height, width))
        self.v = numpy.zeros((height, width))
        self.p = numpy.zeros((height, width)) 

        self.c = c
        self.rho = rho
        self.nu = nu
        self.dt = dt