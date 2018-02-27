'''
Initialization module for CMSE 890 group project.
(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek

History:
v0.1    - LS, 2018-02 -- init
'''

# Import Modules
import numpy as np
import random

def initilization(N_part,spacing):
    '''
    Objective: Initialize the positions and momenta of N particles in the 3D
    
    Inputs: 
            N_part: Number or particles a 1 by 3 array (x,y,z)
            spacing: Particle spacing a 1 by 3 array (x,y,z)
    
    Outputs: 
             position: Position of particles  a N_part by 3 array (x,y,z)
             momenta: Momenta of the particles a N_part by 3 array (x,y,z)
             L: Length of the simulation cell in each direction a 1 by 3 array (x,y,z)
    
    '''
    
    # Number of particles in each direction
    Nx, Ny, Nz = N_part
    
    # Total number of particles
    N_total = Nx*Ny*Nz
    
    # Spacing between partciels in each direction
    sx, sy, sz = spacing

    # Specify the length of the box in each direction
    Lx = sx*Nx
    Ly = sy*Ny
    Lz = sz*Nz

    # Create array of particle position in each direction
    # then shift particles by half of the spacing amount
    # to position particles starting at spacing/2 to L - spacing/2
    x = np.arange(0, Lx, sx) + sx/2
    y = np.arange(0, Ly, sy) + sy/2
    z = np.arange(0, Lz, sz) + sz/2

    # Create a lattice of N_part particles 
    X,Y,Z = np.meshgrid(x,y,z)

    # perfect lattice
    random.seed(1)
    x = X.ravel()
    random.seed(2)
    y = Y.ravel()
    random.seed(3)
    z = Z.ravel()

    if 1 == 1:        
        # Add a perturbation to particles (based on different seeds)
        x += np.random.randn(N_total)*.01
        y += np.random.randn(N_total)*.01
        z += np.random.randn(N_total)*.01
    
    # Store the pariticle positons in array
    position = np.empty([N_total,3])
    position[:,0] = x
    position[:,1] = y    
    position[:,2] = z
    
    # Initialize Momenta from normal distribution 
    mean = 0 # Mean momenta
    width = 1 # Width of distribution
    
    # Store momenta for each particle
    momenta = np.empty([N_total,3])
    momenta[:,0] = np.random.normal(mean,width,N_total)
    momenta[:,1] = np.random.normal(mean,width,N_total)
    momenta[:,2] = np.random.normal(mean,width,N_total)
    
    # Store length of the system in L array
    L = np.array([Lx,Ly,Lz])
    
    return position, momenta, L
