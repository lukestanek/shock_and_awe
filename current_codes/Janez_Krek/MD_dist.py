'''

'''
import numpy as np
import sys
import scipy.constants as co
import copy
import numba
import time
 
class Distribution(object):
   def setvalue(self, N, valx=None, valy=None, valz=None):
      '''
      Initialize vector in 3D of size N with uniform value 'val' or 0 if val is not given.
      
      Args:
         N (int): number of particles in each direction
         valx (various): initial value of [0] cell of each particle
         valy (various): initial value of [1] cell of each particle
         valz (various): initial value of [2] cell of each particle
         
      Returns:
         (array; N x 3): returning vector
      '''
      x = np.zeros((N,3))
      # initialize with given values
      x[:,0] = 0 if valx is None else valx
      x[:,1] = 0 if valy is None else valy
      x[:,2] = 0 if valz is None else valz
      
      return x

   def setvalue_random_direction(self, N, val=0.0):
      '''
      Initialize vector in 3D of size N with uniform value 'val' or 0 if val is not given.
      
      Args:
         N (int): number of particles in each direction
         valx (various): initial value of each particle
         
      Returns:
         (array; N x 3): returning vector
      '''
      x = np.zeros((N,3))
      angle = np.random.random(N)
      
      for i in range(N):
         theta = angle[i]
         phi = angle[-i] 
         x[i] = val*np.array([np.sin(phi)*np.cos(theta), val*np.sin(phi)*np.sin(theta), val*np.cos(phi)])

      return x



   def lattice(self, N, L, margin=0):
      '''
      Initialize vector in 3D of size N with uniform distribution for given length.
      
      Args:
         N (double): total number of particles in all dimensions
         L (double): length in each dimension
         
      Returns:
         (array; N*N*N x 3): returning vector
      '''
      M = int(round(N**(1./3.)))
      x = np.zeros((N,3))
      
      if margin:
         # make effective space for positioning particle smaller forgiven margin
         L  -= 2*margin

      xpos = np.linspace(0, L, M)
      i = 0;
      
      for xi in range(M):
         for yi in range(M):
            for zi in range(M):
               x[i] = [xpos[xi], xpos[yi], xpos[zi]]
               i += 1
      
      return x

   def uniform(self, N, L, minDistance):
      '''
      Initialize vector in 3D of size N with uniform distribution for given length.
      
      Args:
         N (double): number of particles
         L (double): length in each dimension
         minDistance (double): minimum distance between particles
         
      Returns:
         (array; N x 3): returning vector
      '''
      x = np.zeros((N,3))
      x[0] = np.random.uniform(0,L, size=(3))
      i = 1
      
      while i < N:
         p = np.random.uniform(0, L, size=(3))
         tmpMinDist = L
         
         for j in range(i):
            # check distances to all particles
            d = np.sqrt(np.sum((x[j]-p)**2))
            tmpMinDist = min(tmpMinDist, d)
         
         if tmpMinDist > minDistance:
            x[i] = p
            i += 1
      
      return x

   def normal(self, N, sigma=0.01, scale=1):
      '''
      Initialize vector in 3D of size N with normal distribution for given length.
        
      Args:
         N (double): number of particles
         sigma (double): the width of the Normal distribution used to generate values
         scale (double or numpy array): factor for scaling each particle (same for all three dimensions); could be particle mass
         
      Returns:
         (array; N x 3): returning vector
      '''
      x = np.zeros((N,3))

      for i in [0,1,2]:
         # assign velocities in all three coordinates
         x[:,i] = np.random.normal(0, scale=sigma, size=(N))
         x[:,i] = x[:,i]/scale
      
      return x

   def plot(self, x, N=300) :
      '''
      Plot given vector as smoothed distribution (histogram). x has to be 1D numpy vector.
      Method uses 'histogram' and BSpline smoothing to display data gathered form histogram method. 
      Distribution is rendered as curve. 
      
      Args:
         x (array or doubles): values to plot
         N (int): number of points for plotting 
      '''
      import matplotlib.pyplot as plt
   
      plt.figure()
      plt.grid(True) 
    
      bins,  bin_edges = np.histogram(x, bins='auto', density=True)
      num_bins = len(bins)
      xAxis = [0.5*(bin_edges[j]+bin_edges[j+1]) for j in range(num_bins)]
      # smoothed curve
      import scipy.interpolate as interpolate
      from scipy.interpolate import BSpline
      t, c, k = interpolate.splrep(xAxis, bins)
      xnew = np.linspace(np.min(xAxis),np.max(xAxis),N) #300 represents number of points to make between T.min and T.max
      spline = interpolate.BSpline(t, c, k)
      plt.plot(xnew, spline(xnew), linestyle='solid', linewidth=1)
      plt.show()
   
      return 0
