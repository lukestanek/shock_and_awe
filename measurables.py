'''
Module for computing derived quantities (measurables) in the system.

(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek
'''
import numpy as np
from numba import jit

@jit()
def calc_temp(momentum):
    """
    This function calculates the temperature of the system.
    
    Input:
      momentum (numpy array with size num_particles x 3):
            A numpy array containing the momentum of all the particles.
            
    Output:
      Temperature (float):
            A single value that gives the temperature of the system.
    
    """
    Temperature = np.mean(momentum**2)/3
    return Temperature
  

@jit()
def calc_kinetic(momentum):
    """
    This function calculates the kinetic energy of the system.
    
    Input:
      momentum (numpy array with size num_particles x 3):
            A numpy array containing the momentum of all the particles.
            
    Output:
      kinetic (float):
            A single value that gives the kinetic energy of the system.
    
    """
    kinetic = 0
    size = np.size(momentum, axis=0)
    for i in range(0, size):
        # changed "p = np.linalg..." to "v = np.linalg...."
        v = np.linalg.norm(momentum[i]-0)
        kinetic += 0.5*v**2
    return kinetic
  
  
@jit()
def calc_pressure(force, A):
    """
    This function calculates the pressure of the system.
    
    Inputs:
      
      force (numpy array with size num_particles x 3):
            A numpy array containing the force applied all of the particles.
            
      A (float):
            The area of an individual particle.
            
    Output:
      pressure (numpy array with size num_particles x 3):
            A numpy array that contains the pressure each particle is applying.
            
    """
    pressure = force/A
    return pressure
