'''

'''
import numpy as np
import sys
import scipy.constants as co
import copy
import numba
import time

   
#
# === start of definition of various systems (potentials)
#
@numba.jit(nopython=True)
def F(x, L, eps, sigma, rcut, pd, useCutoff=False):
   '''
   Source term;
  
   Args:
      x (numpy array): positions of particles
      m (numpy array): masses of particles
      eps (double): epsilon
      sigma (double): sigma
      rcut (double): cutoff radius
      useCutoff (bool): use cutoff or not
      F (numpy array): returning array for momentum
   '''
   smallDistance = 5e-1
   N = len(x)
   r2cutoff = rcut**2 
   pd[:] = 0
   Lhalf = L/2.

   if useCutoff:
      # use cutoff radius
      for i in xrange(0, N):
         for j in xrange(i+1, N):
            # compute pairwise potentials
            dl = x[i]-x[j]

            for k in xrange(3):
               if dl[k] < Lhalf: dl[k] += Lhalf
               if dl[k] > Lhalf: dl[k] -= Lhalf
            
            l2 = np.sum(dl**2)

            if l2 < r2cutoff:
#                l2 = smallDistance if l2 < smallDistance else l2
               pi = 6*l2**(-4)*(2*l2**(-3) - 1)
               pd[i] += pi*dl
               pd[j] += -pi*dl
            
   else:
      for i in xrange(0, N):
         for j in xrange(i+1, N):
            # compute pairwise potentials
            dl = x[i]-x[j]

            for k in xrange(3):
               if dl[k] < Lhalf: dl[k] += Lhalf
               if dl[k] > Lhalf: dl[k] -= Lhalf
            
            l2 = np.sum(dl**2)
            
#             l2 = smallDistance if l2 < smallDistance else l2
            pi = 6*l2**(-4)*(2*l2**(-3) - 1)
            pd[i] += pi*dl
            pd[j] += -pi*dl
     
   return 


@numba.jit(nopython=True)
def energy(x, p, m, eps, sigma, KE, PE):
   '''
   Compute kinetic, potential and total energies for given positions and moments.

   Args:
      x (double or numpy array): positions
      p (double or numpy array): momentum
      m (numpy array of doubles): particle masses
      eps (double): epsilon
      sigma (double): sigma
      KE (double): kinetic energy - returning array
      PE (double): potential energy - returning array

   '''
   smallDistance = 5e-1
   N = len(x)
   KE[:] = 0
   PE[:] = 0
   
   for i in xrange(N):
      # add together kinetic and potential energy of particles
      KE[i] += np.sum(p[i]**2)
   
      for j in xrange(i+1, N):
         # compute pairwise potentials
         l2 = np.sum((x[i]-x[j])**2)
#          l2 = smallDistance if l2 < smallDistance else l2
         PE[i] += l2**(-6) - l2**(-3)

   KE *= 4*eps/2.
   PE *= 4*eps
   
   return

