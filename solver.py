from os import uname
import numpy

def backwardsDiffX(thing,x,y):
    return thing[y,x]-thing[y,x-1]

def backwardsDiffY(thing,x,y):
    return thing[y,x]-thing[y-1,x]

def sOrderDiffX(thing,x,y):
    return thing[y,x+1]-(2*thing[y,x])+thing[y,x-1]

def sOrderDiffY(thing,x,y):
    return thing[y+1,x]-(2*thing[y,x])+thing[y-1,x]

def centralDiffX(thing,x,y,const):
    return thing[y,x+1]-(thing[y,x-1]*const)

def centralDiffY(thing,x,y,const):
    return thing[y+1,x]-(thing[y-1,x]*const)

def poissonFirst(pMesh,mesh,i,j):
    return ( ((centralDiffX(pMesh,i,j,-1) * mesh.dy**2 + centralDiffY(pMesh,i,j,-1) * mesh.dx**2)) / (2*(mesh.dx**2 + mesh.dy**2)) )

def poissonSecond(mesh,i,j):
    a = mesh.rho#*(mesh.dx**2)*(mesh.dy**2) / (2*((mesh.dx**2)+(mesh.dy**2)))
    b = (((centralDiffX(mesh.u,i,j,1)/(2*mesh.dx)) + (centralDiffY(mesh.v,i,j,1)/(2*mesh.dy)))/mesh.dt)
    c = ((centralDiffX(mesh.u,i,j,1))/2*mesh.dx)**2
    d = 2 * (centralDiffY(mesh.u,i,j,1)/(2*mesh.dy)) * (centralDiffX(mesh.v,i,j,1)/(2*mesh.dx))
    e = (centralDiffY(mesh.v,i,j,1)/(2*mesh.dy)) ** 2

    f = a * (b - c - d - e)

    return f

def solveU(mesh):
    un = mesh.u.copy()

    row, col = mesh.u.shape
    
    for j in range(1,row-1):
            for i in range(1,col-1):
                a = mesh.u[j,i] * backwardsDiffX(mesh.u,i,j) * mesh.dt / mesh.dx 
                b = mesh.v[j,i] * backwardsDiffY(mesh.u,i,j) * mesh.dt / mesh.dy
                c = centralDiffX(mesh.p,i,j,1) * mesh.dt / (mesh.rho * 2 * mesh.dx) 
                d = mesh.nu * (mesh.dt * sOrderDiffX(mesh.u,i,j)/(mesh.dx ** 2) + mesh.dt * sOrderDiffY(mesh.u,i,j)/(mesh.dy ** 2))

                un[j,i] = mesh.u[j,i] - a - b - c + d
    '''
    un[1:-1, 1:-1] = (mesh.u[1:-1, 1:-1] -
                         (mesh.u[1:-1, 1:-1] * mesh.dt / mesh.dx *
                        (mesh.u[1:-1, 1:-1] - mesh.u[1:-1, 0:-2])) -
                         (mesh.v[1:-1, 1:-1] * mesh.dt / mesh.dy *
                        (mesh.u[1:-1, 1:-1] - mesh.u[0:-2, 1:-1])) -
                         (mesh.dt / (2 * mesh.rho * mesh.dx) * (mesh.p[1:-1, 2:] - mesh.p[1:-1, 0:-2])) +
                         (mesh.nu * (mesh.dt / mesh.dx**2 *
                        (mesh.u[1:-1, 2:] - 2 * mesh.u[1:-1, 1:-1] + mesh.u[1:-1, 0:-2]) +
                         mesh.dt / mesh.dy**2 *
                        (mesh.u[2:, 1:-1] - 2 * mesh.u[1:-1, 1:-1] + mesh.u[0:-2, 1:-1]))))
    '''
    return un

def solveV(mesh):
    vn = mesh.v.copy()

    row, col = mesh.v.shape
    
    for j in range(1,row-1):
            for i in range(1,col-1):
                a = mesh.u[j,i] * backwardsDiffX(mesh.v,i,j) * mesh.dt / mesh.dx 
                b = mesh.v[j,i] * backwardsDiffY(mesh.v,i,j) * mesh.dt / mesh.dy
                c = centralDiffY(mesh.p,i,j,1) * mesh.dt / (mesh.rho * 2 * mesh.dy) 
                d = mesh.nu * (mesh.dt * sOrderDiffX(mesh.v,i,j)/(mesh.dx ** 2) + mesh.dt * sOrderDiffY(mesh.v,i,j)/(mesh.dy ** 2))

                vn[j,i] = mesh.v[j,i] - a - b - c + d
    '''
    vn[1:-1,1:-1] = (mesh.v[1:-1, 1:-1] -
                        mesh.u[1:-1, 1:-1] * mesh.dt / mesh.dx *
                       (mesh.v[1:-1, 1:-1] - mesh.v[1:-1, 0:-2]) -
                        mesh.v[1:-1, 1:-1] * mesh.dt / mesh.dy *
                       (mesh.v[1:-1, 1:-1] - mesh.v[0:-2, 1:-1]) -
                        mesh.dt / (2 * mesh.rho * mesh.dy) * (mesh.p[2:, 1:-1] - mesh.p[0:-2, 1:-1]) +
                        mesh.nu * (mesh.dt / mesh.dx**2 *
                       (mesh.v[1:-1, 2:] - 2 * mesh.v[1:-1, 1:-1] + mesh.v[1:-1, 0:-2]) +
                        mesh.dt / mesh.dy**2 *
                       (mesh.v[2:, 1:-1] - 2 * mesh.v[1:-1, 1:-1] + mesh.v[0:-2, 1:-1])))
    '''
    return vn


def solvePressure(mesh,nit):

    b = numpy.empty_like(mesh.p)
    pn = numpy.empty_like(mesh.p)

    row, col = mesh.p.shape
    
    for z in range(1,row-1):
        for s in range(1,col-1):
            b[z,s] = poissonSecond(mesh,s,z)
    '''
    b[1:-1, 1:-1] = (mesh.rho * (1 / mesh.dt * 
                    ((mesh.u[1:-1, 2:] - mesh.u[1:-1, 0:-2]) / 
                     (2 * mesh.dx) + (mesh.v[2:, 1:-1] - mesh.v[0:-2, 1:-1]) / (2 * mesh.dy)) -
                    ((mesh.u[1:-1, 2:] - mesh.u[1:-1, 0:-2]) / (2 * mesh.dx))**2 -
                      2 * ((mesh.u[2:, 1:-1] - mesh.u[0:-2, 1:-1]) / (2 * mesh.dy) *
                           (mesh.v[1:-1, 2:] - mesh.v[1:-1, 0:-2]) / (2 * mesh.dx))-
                          ((mesh.v[2:, 1:-1] - mesh.v[0:-2, 1:-1]) / (2 * mesh.dy))**2))
    '''
    for q in range(nit):
        
        pn = mesh.p.copy()
        
        row, col = mesh.p.shape

        for j in range(1,row-1):
            for i in range(1,col-1):
                mesh.p[j,i] = poissonFirst(pn,mesh,i,j) - mesh.dx**2 * mesh.dy**2 / (2 * (mesh.dx**2 + mesh.dy**2)) * b[j,i]
        '''

        mesh.p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * mesh.dy**2 + 
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * mesh.dx**2) /
                          (2 * (mesh.dx**2 + mesh.dy**2)) -
                          mesh.dx**2 * mesh.dy**2 / (2 * (mesh.dx**2 + mesh.dy**2)) * 
                          b[1:-1,1:-1])
        '''
        mesh.p[:, -1] = mesh.p[:, -2] 
        mesh.p[0, :] = mesh.p[1, :]   
        mesh.p[:, 0] = mesh.p[:, 1]   
        mesh.p[-1, :] = 0 
        




def iterate(mesh,nt,nit):
    for n in range(nt):

        solvePressure(mesh,nit)
        
        un = solveU(mesh)
        vn = solveV(mesh)

        mesh.u = un
        mesh.v = vn

        mesh.u[0, :]  = 0
        mesh.u[:, 0]  = 0
        mesh.u[:, -1] = 0
        mesh.u[-1, :] = 1    
        mesh.v[0, :]  = 0
        mesh.v[-1, :] = 0
        mesh.v[:, 0]  = 0
        mesh.v[:, -1] = 0
        
