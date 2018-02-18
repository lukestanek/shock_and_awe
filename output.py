'''
Module for output - preparing the data into format suitable for visualization or creating visualization.

(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek
'''
import numpy as np
from numba import jit

def write_pos_vel_hist(filename, pos, mom, KE):
   '''
   Write positions and velocities into given file.
   
   Args:
      filename (string): output file name
      pos (numpy array; M x (N,3) ): history of positions 
      mom (numpy array; M x (N,3) ): history of momenta
      mom (numpy array; M x (N,3) ): history of momenta
      KE (numpy array; M x N ): history of momenta
      
   Returns:
      ---
   '''
   Msteps = len(pos)
   N = len(pos[0])

   fp = open(filename, "wb")
   
   print(np.shape(pos[0]), np.shape(mom[0]), np.shape(KE[0]))
   for i in range(Msteps):
      # write current step to file
      np.savetxt(fp, np.c_[pos[i], mom[i], KE[i]], 
                 fmt="%f %f %f %f %f %f %f", header=str(N)+"\n"+ "x y z vx vy vz KE", comments = "")
   
   fp.close()
   print("Positions and momenta saves into file: {0}".format(filename))

   return
      
      
       
