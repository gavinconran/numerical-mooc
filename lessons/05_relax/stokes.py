from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot, cm
from math import pi
import numpy
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

from laplace_helper import L2_rel_error

nx = 41
ny = 41

l = 1.
h = 1.

dx = l/(nx-1)
dy = h/(ny-1)

l1_target = 1e-6
l2_target = 1e-6


# Initialization
omega  = numpy.zeros((ny,nx))
psi  = numpy.zeros((ny,nx))
U_top  = numpy.ones((ny-3))
x = numpy.linspace(0,1,nx)
y = numpy.linspace(0,1,ny)
psi[-2,:] = 1


def L1norm(new, old):
    norm = numpy.sum(numpy.abs(new-old))
    return norm

def L2_error(p, pn):
    return numpy.sqrt(numpy.sum((p - pn)**2)/numpy.sum(pn**2))

def laplace2d(p, psi, l2_target):
    '''Iteratively solves the Laplace equation using the Jacobi method
    
    Parameters:
    ----------
    p: 2D array of float
        Initial potential distribution
    l2_target: float
        target for the difference between consecutive solutions
        
    Returns:
    -------
    p: 2D array of float
        Potential distribution after relaxation
    '''
    
    l2norm = 1
    pn = numpy.empty_like(p)
    iterations = 0
    while l2norm > l2_target:
        pn = p.copy()
        p[1:-1,1:-1] = .25 * (pn[1:-1,2:] + pn[1:-1, :-2] \
                              + pn[2:, 1:-1] + pn[:-2, 1:-1])

        # update velocity
        U_top = (psi[-2,:] - psi[-1,:]) /(dy)
        U_bottom = (psi[1,:] - psi[0,:]) /(dy)
        U_left = (psi[:,-2] - psi[:,-1]) /(dy)
        U_right = (psi[:,1] - psi[:,0]) /(dy)

        
        # Apply Neumann BC to omega
        p[-1,2:-1] = (1/2*dy**2) * (8*psi[-1,1:-2] - psi[-1,:-3]) #- 3*(U_top[2:-1]) / dy
        p[0,2:-1] = (1/2*dy**2) * (8*psi[0,1:-2] - psi[0,:-3]) #- 3*U_bottom[2:-1] / dy
        #p[2:-1,0] = (1/2*dx**2) * (8*psi[1:-2,0] - psi[:-3,0]) #- 3*U_left[2:-1] / dx
        #p[2:-1,-1] = (1/2*dx**2) * (8*psi[1:-2,-1] - psi[:-3,-1]) #- 3*U_right[2:-1] / dx
        #p[1:-1, -1] = p[1:-1, -2]
        l2norm = L2_error(p, pn)
     
    return p

def poisson_2d(omega, psi, dx, dy, l1_target):
    '''Performs Jacobi relaxation
    
    Parameters:
    ----------
    omega : 2D array of floats
        Vorticity
    psi : 2D array of floats
        Stream function
    dx: float
        Mesh spacing in x direction
    dy: float
        Mesh spacing in y direction
    l2_target: float
        Target difference between two consecutive iterates
    
    Returns:
    -------
    omega: 2D array of float
        Distribution after relaxation
    psi: 2D array of float
        Distribution after relaxation    
    '''

    l1_norm = 1
    iterations = 0
    l1_conv = []
    
    omega_old = omega.copy()
    
    while l1_norm > l1_target:
        
        # compute Psi
        psi_d = psi.copy()
        
        psi[1:-1,1:-1] = 1/(2*(dx**2 + dy**2)) * \
                        ((psi_d[1:-1,2:] + psi_d[1:-1,:-2])*dy**2 +\
                        (psi_d[2:,1:-1] + psi_d[:-2,1:-1])*dx**2 -\
                         omega[1:-1,1:-1])
        
        
        # Dirichlet BCs
        psi[0,:] = 0    #psi[1,:]
        psi[-1,:] = 0   #psi[-2,:] + 1
        psi[:,0] = 0    #psi[:,1]
        psi[:,-1] = 0   #psi[:,-2]
        
        
        # computr Omega
        omega = laplace2d(omega.copy(), psi.copy(), l2_target)
                

        
        # complute L1 norm
        l1_norm = L1norm(psi,psi_d)
        
        iterations += 1
        l1_conv.append(l1_norm)
    
    print('Number of Jacobi iterations: {0:d}'.format(iterations))
    return omega, psi, l1_conv        

# compute psi
#print(psi) 
omega, psi, l2_conv = poisson_2d(psi.copy(), omega.copy(), dx, dy, l1_target)


#pyplot.figure(figsize=(7,7))
#pyplot.contourf(x,y,omega,20, cmap=cm.viridis);

pyplot.figure(figsize=(7,7))
pyplot.contourf(x,y,psi,20, cmap=cm.viridis);
pyplot.show()
print('max omega: ', numpy.amax(numpy.abs(omega)))
print('max psi: ', numpy.amax(numpy.abs(psi)))

print(numpy.round(psi[32,::8], 4)) 
#print(psi) 


