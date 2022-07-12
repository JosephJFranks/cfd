import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

from setup import Mesh
from solver import iterate

nx = 41 
ny = 41
nt = 500
nit = 50
c = 1

rho = 1
nu = .1
dt = .001

Fluid = Mesh(nx,ny,c,nu,rho,dt)

#self,width,height,c,nu,rho,dt

iterate(Fluid,nt,nit)

fig = pyplot.figure(figsize=(11,7), dpi=100)
# plotting the pressure field as a contour
pyplot.contourf(Fluid.X, Fluid.Y, Fluid.p, alpha=0.5, cmap=cm.viridis)  
pyplot.colorbar()
# plotting the pressure field outlines
pyplot.contour(Fluid.X, Fluid.Y, Fluid.p, cmap=cm.viridis)  
# plotting velocity field
pyplot.quiver(Fluid.X[::2, ::2], Fluid.Y[::2, ::2], Fluid.u[::2, ::2], Fluid.v[::2, ::2]) 
pyplot.xlabel('X')
pyplot.ylabel('Y');

pyplot.show()
