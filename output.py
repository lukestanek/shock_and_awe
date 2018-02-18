'''
Module for output - preparing the data into format suitable for visualization or creating visualization.

(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek
'''
import numpy as np
from numba import jit

def write_pos_vel_hist(filename, pos, mom, KE, pis, Lx, Ly, Lz):
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
   tpos = np.zeros((8,3))
   p_tags = [0]*N
   b_tags = [1]*8
   fp = open(filename, "wb")
   
   print(np.shape(pos[0]), np.shape(mom[0]), np.shape(KE[0]))
   for i in range(Msteps):
	
      # pistion postion
      tpos[0] = [0,0,pis[i]]
      tpos[1] = [0,Ly,pis[i]]
      tpos[2] = [Lx,Ly,pis[i]]
      tpos[3] = [Lx,0,pis[i]]
	  
      # right boundary
      tpos[4] = [0,0,Lz]
      tpos[5] = [0,Ly,Lz]
      tpos[6] = [Lx,Ly,Lz]
      tpos[7] = [Lx,0,Lz]
      # write current step to file
      np.savetxt(fp, np.c_[p_tags, pos[i], mom[i], KE[i]], 
                 fmt="%f %f %f %f %f %f %f %f", header=str(N+8)+"\n"+ "type x y z vx vy vz KE", comments = "")
   
      np.savetxt(fp, np.c_[b_tags, tpos, np.zeros((8,3)), np.zeros(8)], 
                 fmt="%f %f %f %f %f %f %f %f", header="", comments = "")
   fp.close()
   print("Positions and momenta saves into file: {0}".format(filename))

   return
      
      
       
