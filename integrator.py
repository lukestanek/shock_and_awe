'''
Module for integrator.
  
(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek

History:
v0.1    - DL, 2018-02 -- init
v0.2    - JK, 2018-02-22 -- added (nopython=True) to jit()
v0.3    - JK, 2018-02-22 -- small rewrite in calc_force to make it faster + using prange
'''  
import numpy as np
from numba import jit, prange
import boundaries
import sys

@jit()#nopython=True)#, parallel=True)
def vel_ver(position, momentum, Piston_p, Piston_Momentum, dt, force, x_len, y_len, Mirror_Position, radius):
    """
    This is the function that will update our position and momentum arrays. It
    assumes a rectangular piston cross-section and it applies boundary conditions
    using the boundary.py functions. The force function it calls is defined below
    it. The force that the function outputs will be the input force for the next
    time step.
    
    Inputs:
      position (numpy array with size num_particles x 3):
            A numpy array containing the positions of all the particles
            
      momentum (numpy array with size num_particles x 3):
            A numpy array containing the momentum of all the particles
            
      Piston_p (float):
            The current position of the piston
            
      Piston_Momentum (float):
            Momentus of the piston
            
      dt (float):
            The size of the timestep
      
      force (numpy array with size num_particles x 3):
            A numpy array containing the force vector for each particle
            
      x_len (float):
            The length of the box in the x direction. Assumes that the
            box starts at x = 0 and goes in the positive direction
    
      y_len (float):
             The length of the box in the y direction. Assumes that the
             box starts at y = 0 and goes in the positive direction
          
      Mirror_Position (float):
            Position of the back of the cell. Assumes that the box starts
            at the z = 0 and goes in the positive direction
            
      radius (float):
            The radius of the sphere that each particle interacts with
            
    Outputs:
      position (numpy array with size num_particles x 3):
             A numpy array containing the new positions of all the particles
             
      momentum (numpy array with size num_particles x 3):
             A numpy array containing the new momentum of all the particles
             
      force (numpy array with size num_particles x 3):
             A numpy array containing the new force vector for each particle
             
    """
    # Updates the halfstep momentum and the unedited position
    momentum += 0.5*dt*force
    position += momentum*dt

    # Applies the boundary conditions to the x, y positions
    #size = np.size(position, axis=0)
    size = len(position)        # JK, 2018-02-22, nopython=True
    position = boundaries.periodic_boundary_position(position , size, x_len, y_len)
    
    # Applies the momentum mirror to the z positions
    position, momentum = boundaries.Momentum_Mirror(position, momentum, Piston_Momentum, Piston_p, Mirror_Position, dt, size)

    # Updates the final force and momentum
    force = calc_force(position, radius, x_len, y_len)
    momentum += 0.5*dt*force

    return position, momentum, force
  
  
@jit()#nopython=True)#, parallel=True)
def calc_force(position, radius, x_len, y_len):
    """
    This is the function that will calculate the forces for the 
    velocity verlay integrator above. Like the integrator, it assumes
    a rectangular cross section and applies those boundary conditions
    using the functions defined in boundary.py.
    
    Inputs:
      position (numpy array with size num_particles x 3):
            A numpy array containing the positions of all the particles
            
      radius (float):
            The radius of the sphere of interaction for the force calculation
                    
      x_len (float):
            The length of the box in the x direction. Assumes that the
            box starts at x = 0 and goes in the positive direction.
    
      y_len (float):
             The length of the box in the y direction. Assumes that the
             box starts at y = 0 and goes in the positive direction.
            
    Outputs:
      force (numpy array with size num_particles x 3):
            A numpy array containing the force vector for each particle
    """
    
    # Creates force array
    size = len(position)        # JK, 2018-02-22; because of noptyhon=True
    force = np.zeros((size, 3))
    radius2 = radius**2

    # Calculates distances arrays using the periodic boundary conditions
    x_diff, y_diff, z_diff = boundaries.periodic_boundary_force(position, size, x_len, y_len)
    r_tilde = x_diff**2 + y_diff**2 + z_diff**2

    for i in prange(0,size-1):
        for j in prange(i+1,size):
            # using prange for paralllel 
            # If particle is in poor man's radius, calculates force
            if r_tilde[i][j] <= radius2:
                S = 6*( 2*(r_tilde[i][j]**-7) - (r_tilde[i][j]**-4) )
                Sx,Sy,Sz = S*x_diff[i][j],S*y_diff[i][j],S*z_diff[i][j]

                force[i][0] += Sx
                force[j][0] -= Sx

                force[i][1] += Sy
                force[j][1] -= Sy

                force[i][2] += Sz
                force[j][2] -= Sz

    return force

@jit()#nopython=True)#, parallel=True)
def calc_potential_energy(position, radius, x_len, y_len):
    """
    This is the function that will calculate total potential energy in the
    system.
    
    Inputs:
      position (numpy array with size num_particles x 3):
            A numpy array containing the positions of all the particles
            
      radius (float):
            The radius of the sphere of interaction for the force calculation
                    
      x_len (float):
            The length of the box in the x direction. Assumes that the
            box starts at x = 0 and goes in the positive direction.
    
      y_len (float):
             The length of the box in the y direction. Assumes that the
             box starts at y = 0 and goes in the positive direction.
            
    Outputs:
      PE (float):
            Total potential energy in the system
    """
    # Creates force array
    size = len(position)        # JK, 2018-02-22; because of noptyhon=True
    radius2 = radius**2
    PE = 0.0

    # Calculates distances arrays using the periodic boundary conditions
    x_diff, y_diff, z_diff = boundaries.periodic_boundary_force(position, size, x_len, y_len)
    r_tilde = x_diff**2 + y_diff**2 + z_diff**2

    for i in range(0,size-1):
        for j in range(i+1,size):
            # If particle is in poor man's radius, calculates force
            if r_tilde[i][j] <= radius2:
                PE += r_tilde[i][j]**(-6) - r_tilde[i][j]**(-3)

    return PE

