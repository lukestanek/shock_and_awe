'''
Module for integrator.
  
(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek
'''  
import numpy as np
from numba import jit

@jit()
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
            
      dt (float):
            The size of the timestep
      
      force (numpy array with size num_particles x 3):
            A numpy array containing the force vector for each particle
            
      x_len (float):
            The length of the box in the x direction. Assumes that the
            box starts at x = 0 and goes in the positive direction.
    
      y_len (float):
             The length of the box in the y direction. Assumes that the
             box starts at y = 0 and goes in the positive direction.
            
    Outputs:
      new_position (numpy array with size num_particles x 3):
              A numpy array containing the new positions of all the particles
      new_momentum (numpy array with size num_particles x 3):
              A numpy array containing the new momentum of all the particles
      new_force (numpy array with size num_particles x 3):
              A numpy array containing the new force vector for each particle
    """
    # Updates the Piston Position
    Piston_p = Piston_Position(Piston_p, Piston_Momentum, dt)
    
    # Updates the halfstep momentum and the unedited position
    momentum_half = 0.5*dt*force + momentum
    new_position = momentum_half*dt + position

    # Applies the boundary conditions to the x, y positions
    size = np.size(position, axis=0)
    new_position = periodic_boundary_position(new_position , size, x_len, y_len)
    
    # Applies the momentum mirror to the z positions
    new_position = Momentum_Mirror(new_position, momentum_half, Piston_Momentum, Piston_Position, Mirror_Position, size)

    # Updates the final force and momentum
    new_force = calc_force(new_position, radius, x_len, y_len)
    new_momentum = (0.5*dt*new_force + momentum_half)

    return new_position, new_momentum, new_force
  
  
@jit()
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
    size = np.size(position, axis=0)
    force = np.zeros((size, 3))

    # Calculates distances arrays using the periodic boundary conditions
    x_diff, y_diff, z_diff = periodic_boundary_force(position, size, x_len, y_len)
    r_tilde = x_diff**2 + y_diff**2 + z_diff**2
    
    for i in range(size-1):
        for j in range(i+1,size):
          
            # If particle is in poor man's radius, calculates force
            if r_tilde[i][j] <= radius**2:
                S = 6*( 2*(r_tilde[i][j]**-7) - (r_tilde[i][j]**-4) )
                Sx,Sy,Sz = S*x_diff[i][j],S*y_diff[i][j],S*z_diff[i][j]

                force[i][0] += Sx
                force[j][0] -= Sx

                force[i][1] += Sy
                force[j][1] -= Sy

                force[i][2] += Sz
                force[j][2] -= Sz

            else:
                force[i][0] += 0
                force[j][0] -= 0

                force[i][1] += 0
                force[j][1] -= 0

                force[i][2] += 0
                force[j][2] -= 0

    return force

