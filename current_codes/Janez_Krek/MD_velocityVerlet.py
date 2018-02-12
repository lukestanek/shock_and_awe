'''

'''
import numpy as np
import sys
import scipy.constants as co
import copy
import numba
import time
from MD_plotting import Plot
from MD_dist import Distribution
import MD_LennardJones as LJ
 
@numba.jit(nopython=True)
def update(dt, Ft, L, eps, sigma, rcutoff, x, v, useCutoff=False):
   '''
   Update positions and velocities using velocty Verlet. Positions and velocities are updated in-line - new values replace old ones.
   '''
   # first part of the push
   v += 0.5*Ft*dt
   x += v*dt
   LJ.F(x, L, eps, sigma, rcutoff, Ft, useCutoff)
   # second part of the push
   v += 0.5*Ft*dt
   return


#@numba.jit('double(intp,intp,double[:,:])',nopython=True)
@numba.jit(nopython=True)
def square_vel_mass_sum(N,m,v):
   #tmpvel = vel**2
   out = 0.0
   for i in xrange(N):
      v2 = 0.0
   
      for j in xrange(3):
         v2 += v[i][j]**2
      
      out += m[i]*v2

   return out

@numba.jit(nopython=True)
def compute_temp(N, m, v):
   # compute <m v^2>
   en = square_vel_mass_sum(N,m,v)
   # temp is 1/3 <m v^2>/N
   T = en/N/3.

   return T


@numba.jit(nopython=True)
def scale(Ttarget, Tcurr, Ncalib, v):
   '''
   Scale velocities with current Temp to target temperature.
   '''
   s = np.sqrt(1 + 1./Ncalib*(Ttarget*1.0/Tcurr-1))
   v *= s
   return
