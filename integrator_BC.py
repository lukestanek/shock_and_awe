'''
Module for integrator.
'''
import numpy as np

@jit()
def vel_ver(position, momentum, dt, force, radius, length):
    """
    This is the function that will update our position and momentum arrays.
    It assumes a square piston cross-section and it includes the boundary
    conditions. The force function it calls is defined below it. The force
    that the function outputs will be the input force for the next time step.

    Inputs:

      position (numpy array with size num_particles x 3):
            A numpy array containing the positions of all the particles

      momentum (numpy array with size num_particles x 3):
            A numpy array containing the momentum of all the particles

      dt (float): The size of the timestep

      force (numpy array with size num_particles x 3):
            A numpy array containing the force vector for each particle

      radius (float): The radius of the sphere of interaction for the 
                    force calculation
                    
      length (float): The side length of the square that makes up the
                    cross section of the piston

    Outputs:

      new_position (numpy array with size num_particles x 3):
            A numpy array containing the new positions of all the particles

      new_momentum (numpy array with size num_particles x 3):
            A numpy array containing the new momentum of all the particles

      new_force (numpy array with size num_particles x 3):
            A numpy array containing the new force vector for each particle
    """
    
    # Updates the halfstep momentum and the unedited position
    momentum_half = 0.5*dt*force + momentum
    new_position = momentum_half*dt/mass + position

    # Applies the boundary conditions to the x, y positions
    new_position[:,0:2] = new_position[:,0:2] % length

    # Updates the final force and momentum
    new_force = calc_force(new_position, radius, length)
    new_momentum = (0.5*dt*new_force + momentum_half)

    return new_position, new_momentum, new_force
  
  
@jit()
def calc_force(position, radius, length):
    """
    This is the function that will calculate the forces for the 
    velocity verlay integrator above. Like the integrator, it assumes
    a square cross section and applies those boundary conditions.

    Inputs:

      position (numpy array with size num_particles x 3):
            A numpy array containing the positions of all the particles

      radius (float): The radius of the sphere of interaction for the 
                    force calculation
                    
      length (float): The side length of the square that makes up the
                    cross section of the piston

    Outputs:

      force (numpy array with size num_particles x 3):
            A numpy array containing the force vector for each particle
    """
    # Creates force array
    size = np.size(position, axis=0)
    force = np.zeros((size, 3))

    for i in range(size-1):
        for j in range(i+1,size):

            # Calculates distances and applies boundary conditions
            x_diff = position[i][0] - position[j][0]
            y_diff = position[i][1] - position[j][1]
            z_diff = position[i][2] - position[j][2]

            x_diff = x_diff - np.floor((x_diff+0.5*length)/length)*length
            y_diff = y_diff - np.floor((y_diff+0.5*length)/length)*length

            r_tilde = x_diff**2 + y_diff**2 + z_diff**2

            # If particle is in poor man's radius, calculates force
            if r_tilde**2 <= radius**2:
                S = 24*epsilon*( 2*(sigma**12)*(r_tilde**-7) - (sigma**6)*(r_tilde**-4) )
                Sx,Sy,Sz = S*x_diff,S*y_diff,S*z_diff

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
'''  
(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek
'''  
