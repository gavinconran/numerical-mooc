import numpy


def f(U,gamma):
    """Returns the right-hand side of the rocket system of equations.
    
    Parameters
    ----------
    u : array of float
        array containing the solution at time n.
        
    Returns
    -------
    F : array of float
        array containing the RHS given u.
    """
    
    # Extract parameters from vector U
    u1 = U[0] # density
    u2 = U[1] # density * velocity
    u3 = U[2] # density * energyT
    u4 = U[3] # presure
    
    # cmpute pressure and energyT
    p = (gamma - 1) * (u3 - 0.5 * u2**2 / u1)
    e = u4 / ((gamma -1)*u1)
    eT = e + 0.5*u1**2
    
    f_vector = numpy.array([ u2,
                       u2**2/u1 + p,
                       (u3 + p) * u2/u1,
                       p,
                       eT])  
    
    return f_vector

# Setup conditions
# Construct grid
nx = 81
dx = 0.25
dt = 0.0002
T = 0.001 #0.01
nt = int(T/dt)
gamma = 1.4
x_L = -10
x_R = 10

# Initial conditions
#Boundary conditions LEFT
rho_L = 1.
vel_L = 0.
pressure_L = 100. * 1000

#Boundary conditions RIGHT
rho_R = 0.125
vel_R = 0.
pressure_R = 10. * 1000

x=numpy.linspace(x_L,x_R,nx)

# Density
rho0=numpy.zeros(nx)
rho0[0:int((nx-1)/2)]=rho_L
rho0[int((nx-1)/2):nx]=rho_R

# Velocity
vel0=numpy.zeros(nx)
vel0[0:int((nx-1)/2)]=vel_L
vel0[int((nx-1)/2):nx]=vel_R

# pressure
pressure0=numpy.zeros(nx)
pressure0[0:int((nx-1)/2)]=pressure_L
pressure0[int((nx-1)/2):nx]=pressure_R

# energy
energy = pressure0 / ((gamma-1)*rho0)
energyT0=numpy.zeros(nx)
energyT0[0:int((nx-1)/2)]=(energy[0:int((nx-1)/2)] + 0.5*vel_L)
energyT0[int((nx-1)/2):nx]=(energy[int((nx-1)/2):nx] + 0.5*vel_R)

## Test function f
U = numpy.zeros((5, len(rho0)))
U[0] = rho0
U[1] = rho0 * vel0
U[2] = rho0 * energyT0
U[3] = pressure0
U[4] = energyT0

