'''
Module for boundaries.

(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek

History:
v0.1    - TD, RM, 2018-02 -- init
v0.2    - JK, 2018-02-22 -- added piston movement end time; added (nopython=True, parallel=True) to jit()
v0.3    - JK, 2018-02-22 -- small rewrite in periodic_boundary_force to make it faster + using prange
'''
import numpy as np
from numba import jit, prange

@jit()#nopython=True)#, parallel=True)
def periodic_boundary_force(pos, n_par, x_len, y_len):
    """
    This function applies the periodic boundaries to the
    particle postions in order to determine the distance
    between all the particles.
    
    Inputs:
    
    pos (2D numpy array size: number of particles x 3): 
    A numpy array that contains the positions of all the
    particles in 3D space.
    
    n_par (float): The number of particles.
    
    x_len (float): The length of the box in the x direction.
                    Assumes that the box starts at x = 0 and goes
                    in the positive direction.
    
    y_len (float): The length of the box in the y direction.
                    Assumes that the box starts at y = 0 and goes
                    in the positive direction.
                    
    Outputs:
    
    x_diff (1D numpy array size: number of particles x number of particles): 
            The difference in x positions between all the particles
            after periodic boundary conditions are applied.
            
    y_diff (1D numpy array size: number of particles x number of particles): 
            The difference in y positions between all the particles
            after periodic boundary conditions are applied.
            
    z_diff (1D numpy array size: number of particles x number of particles): 
            The difference in z positions between all the particles
    """
    
    # Initialize the distances in positions arrays.
    minDistance = 0.6               # JK, 2018-02-24
    # JK, 2018-02-22; nopython=True
    x_diff = np.zeros((n_par,n_par))
    y_diff = np.zeros((n_par,n_par))
    z_diff = np.zeros((n_par,n_par))

    for par_1 in range(n_par):
        for par_2 in range(par_1+1, n_par):
            # loop over all particles; using prange for parallel; JK, 2018-02-22
            # Find x, y, and z distances
            x_diff[par_1, par_2] = pos[par_1, 0] - pos[par_2, 0]
            y_diff[par_1, par_2] = pos[par_1, 1] - pos[par_2, 1]
            z_diff[par_1, par_2] = pos[par_1, 2] - pos[par_2, 2]
            
            # Applies minimum distance boundary
            if abs(x_diff[par_1, par_2]) < minDistance:
                if x_diff[par_1, par_2] > 0:
                    x_diff[par_1, par_2] = minDistance
                else:
                    x_diff[par_1, par_2] = -minDistance
            		
            if abs(y_diff[par_1, par_2]) < minDistance:
                if y_diff[par_1, par_2] > 0:
                    y_diff[par_1, par_2] = minDistance
                else:
                    y_diff[par_1, par_2] = -minDistance
            		
            if abs(z_diff[par_1, par_2]) < minDistance:
                if z_diff[par_1, par_2] > 0:
                    z_diff[par_1, par_2] = minDistance
                else:
                    z_diff[par_1, par_2] = -minDistance
            		
            # Applies periodic boundary condition
            if abs(x_diff[par_1, par_2]) > x_len:
                if x_diff[par_1, par_2] > 0:
                    x_diff[par_1, par_2] -= x_len
                else:
                    x_diff[par_1, par_2] += x_len 
            		
            if abs(y_diff[par_1, par_2]) > y_len:
                if y_diff[par_1, par_2] > 0:
                    y_diff[par_1, par_2] -= y_len
                else:
                    y_diff[par_1, par_2] += y_len 

            # store same distance to "other" particle
            x_diff[par_2, par_1] = -x_diff[par_1, par_2]
            y_diff[par_2, par_1] = -y_diff[par_1, par_2]
            z_diff[par_2, par_1] = -z_diff[par_1, par_2]

    return(x_diff, y_diff, z_diff)

@jit()#nopython=True)#, parallel=True)  
def periodic_boundary_position(pos , n_par, x_len, y_len):
    '''
    This function moves the particles according to the periodic boundaries
    in the x and y directions.
    
    Inputs:
    
    pos (2D numpy array size: number of particles x 3): 
    A numpy array that contains the positions of all the
    particles in 3D space.
    
    n_par (float): The number of particles.
    
    x_len (float): The length of the box in the x direction.
                    Assumes that the box starts at x = 0 and goes
                    in the positive direction.
    
    y_len (float): The length of the box in the y direction.
                    Assumes that the box starts at y = 0 and goes
                    in the positive direction.
                    
    Outputs:
    
    new_pos (2D numpy array size: number of particles x 3): 
    A numpy array that contains the positions of all the
    particles in 3D space after applying periodic boundary
    conditions.
    '''
    # Creates an array to store all the new positions
    # Apply periodic boundary conditions. We are not applying periodic
    # boundaries to the z component
    
    for par in range(n_par):
        # using prange for parallel; JK, 2018-02-22
        # x component
        if pos[par,0] > x_len or pos[par, 0] < 0:
            pos[par,0] = pos[par, 0] % x_len
            
        # y component    
        if pos[par,1] > y_len or pos[par, 1] < 0:
            pos[par,1] = pos[par, 1] % y_len
            
    return pos


#@jit()#nopython=True)#, parallel=True)
def Momentum_Mirror(Position, Momentum, Piston_Momentum, Piston_Position, Mirror_Position, dt, N_Par):
    #This is the function for the momentum mirror and the Piston. It takes in: 
    
    #* Particle_Position = This should be a np.array with dimensions (N, 3).
    #* Particle_Momentum = This should be a np.array with dimensions (N, 3).
    #* Piston_Momentum = This should be a float.
    #* Piston_Position = This will be a float that is updated after every loop.
    #* Mirror_Position = This should also be a float.
    #* dt = time step
    #* Number of particles = This should be a Constant Integer.
    
   # This function will check each particles z position, since that is the dimension we chose, 
   # and update the new z position if the position violates our conditions. Finally the function 
   # returns the Particles Positions and Momentums.
    
    for particle in range(N_Par):
        # using prange for parallel; JK, 2018-02-22
        Particle_Position = Position[particle]
        Particle_Momentum = Momentum[particle]
        
        if Particle_Position[2] < Piston_Position:
            Particle_Position[2] = Piston_Position + (Piston_Position - Particle_Position[2])
            Particle_Momentum[2] = -Particle_Momentum[2] + Piston_Momentum
            
        elif Particle_Position[2] > Mirror_Position:
            Particle_Position[2] = Mirror_Position - (Particle_Position[2] - Mirror_Position)
            Particle_Momentum[2] = -(Particle_Momentum[2])  
            
    return Position, Momentum


@jit()#nopython=True)
def calc_Piston_Position(Piston_Position, Piston_Velocity, dt, pistonEndTime, t):
    #This function is extremely straight forward, the function takes in:
    
    #* Piston_Position = This should be Float.
    #* Piston Velocity = This should be a Float.
    #* dt = This should be Float.
    #* pistonEndTime = end time of piston moving (JK, 2018-02-22)
    #* t = current time in the simulation (JK, 2018-02-22)
    
    #The function takes the piston's position and updates it. Yup.
    
    Piston_Position += Piston_Velocity * dt
    
    return Piston_Position


