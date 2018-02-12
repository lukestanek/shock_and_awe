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
import MD_velocityVerlet as velocity_verlet
import MD_boundaries as boundaries



def runMD(in_x0, in_p0, m, params, usePBC=False, runEQ=False, useCutoff=False, outFile=None, outFileTimeEvol=None, skipSteps=0, progress=True) :
   '''
   Run case with given dt and initial angle given by 'x0'.

   Args:
      in_x0 (array of doubles): initial positions; including walls (data [0] and [-1])
      in_p0 (array of doubles): initial moments; including walls (data [0] and [-1])
      m (numpy array of doubles): particle masses
      params (dict): simulation parameters
      usePBC (bool): use periodic boundaries or not
      runEQ (bool): run equilibrium steps (True) or not (False)
      outFile (string or None): write positions, velocities, energies to given file, XYZ style; all time steps are written to one file; in case of None, no writing
      skipStep (int): number of steps to skip between writes into file
      outFileEnergies (string or None): write energies in every time step to given file
   '''
   maxTime = params['T']
   dt = params['dt']
   L = params['L']
   eps = params['eps']
   sigma = params['sigma']
   Ttarget = params['Ttarget']
   Ncalib = params['Ncalibation']
   rcutoff = params['r_cutoff']
   
   # scale factors
   scaleX = sigma
   scaleE = sigma/2.*np.sqrt(m/eps)
   scaleP = 2*np.sqrt(eps)/m

   fp = None
   fpTime = None
   
   if outFile is not None:
      fp = open(outFile, "w")
   if outFileTimeEvol is not None:
      fpTime = open(outFileTimeEvol, "w")
         
   N = len(in_x0)
   M = int(round(maxTime*1.0/dt))
#    print percent
#    print "Running case for dt={0}, T={1}, M={2}.".format(dt, T, M)
   
   # init arrays for data
   pos = np.zeros((M,N,3))
   vel = np.zeros_like(pos)

   x =  np.zeros((N,3))
   p = copy.deepcopy(x)
   Ft = copy.deepcopy(x)
   KE = np.zeros(N)
   PE = np.zeros(N)
   totalKE = np.zeros(M)
   totalPE = np.zeros(M)
   totalH = np.zeros(M)
   T = np.zeros(M)
   # init position and moment   
   x = in_x0
   p = in_p0

   # initial energy
   LJ.energy(x, p, m, eps, sigma, KE, PE)
   Temp = velocity_verlet.compute_temp(N, m, p)
   pos[0] = x
   vel[0][:,0] = p[:,0]*scaleP
   vel[0][:,1] = p[:,1]*scaleP
   vel[0][:,2] = p[:,2]*scaleP
   totalKE[0] = np.sum(KE*scaleE)
   totalPE[0] = np.sum(PE*scaleE)
   totalH[0] = totalKE[0]+totalPE[0]
   T[0] = Temp
   
#    if fp is not None:
#       # write current step to file
#       np.savetxt(fp, np.c_[pos[0], vel[0], m, KE*scaleE, PE*scaleE, (KE+PE)*scaleE], 
#                  fmt="%f %f %f %f %f %f %f %f %f %f", header=str(N)+"\n"+ "x y z vx vy vz m KE PE H", comments = "")
#       
#    timeF = 0
#    timeE = 0

   if runEQ:
      if progress:
         print "==> Equilibration run...."
         print " - computing...",

      KE[:] = 0
      PE[:] = 0
      totalKE[:] = 0
      totalPE[:] = 0
      totalH[:] = 0
      T[:] = 0
      Meq = M
      percent = [int(el) for el in Meq*np.linspace(0,1.,21)]
      percent.pop(0)
      # compute forces on particles in intial position
      LJ.F(x, L, eps, sigma, rcutoff, Ft, useCutoff=useCutoff)
      
      for i in xrange(1, Meq):
         t = dt*i
         
         if progress:
            if( i == percent[0] ):
               print "{0:d}%...".format(int((i)*100.0/M)),
               percent.pop(0)
               sys.stdout.flush()

         time1 = velocity_verlet.update(dt, Ft, L, eps, sigma, rcutoff, x, p, useCutoff=useCutoff)
         
         if usePBC:
            boundaries.pbc(L, x, p)     # do the PBC
            
         # scale velocities
         Temp = velocity_verlet.compute_temp(N, m, p)
         velocity_verlet.scale(Ttarget, Temp, Ncalib, p)
         LJ.energy(x, p, m, eps, sigma, KE, PE)
   
         pos[i] = x
         vel[i][:,0] = p[:,0]*scaleP
         vel[i][:,1] = p[:,1]*scaleP
         vel[i][:,2] = p[:,2]*scaleP
         totalKE[i] = np.sum(KE*scaleE)
         totalPE[i] = np.sum(PE*scaleE)
         T[i] = Temp

      print "100%. Done!"

      if outFileTimeEvol is not None:
         # save values into file ending with "-eq"
         fn = outFileTimeEvol+"-eq"
         fpTime = open(fn, "w")
      
         fpTime.write("# {0} {1}".format(maxTime, dt))
         np.savetxt(fpTime, np.c_[np.linspace(0,maxTime,M), totalPE[:Meq], totalKE[:Meq], totalH[:Meq], T[:Meq]], fmt="%f %f %f %f %f", header="t PE KE H T")
         fpTime.close()
         print " - energies/temperatures written into: ", fn
   
   else:
      print "Equilibrium steps skipped."


   # init arrays for energies      
   KE[:] = 0
   PE[:] = 0
   totalKE[:] = 0
   totalPE[:] = 0
   totalH[:] = 0
   T[:] = 0
   # store current state as init state before real run
   LJ.energy(x, p, m, eps, sigma, KE, PE)
   Temp = velocity_verlet.compute_temp(N, m, p)
   pos[0] = x
   vel[0][:,0] = p[:,0]*scaleP
   vel[0][:,1] = p[:,1]*scaleP
   vel[0][:,2] = p[:,2]*scaleP
   totalKE[0] = np.sum(KE*scaleE)
   totalPE[0] = np.sum(PE*scaleE)
   totalH[0] = totalKE[0]+totalPE[0]
   T[0] = Temp
   timeF = 0
   timeE = 0
   
   percent = [int(el) for el in M*np.linspace(0,1.,21)]
   percent.pop(0)

   if progress:
      print
      print "==> Simulation run"
      print " - computing...",
      
   # compute forces on particles in intial position
   LJ.F(x, L, eps, sigma, rcutoff, Ft, useCutoff=useCutoff)
   
   for i in xrange(1, M):
      t = dt*i
      
      if progress:
         if( i == percent[0] ):
            print "{0:d}%...".format(int((i)*100.0/M)),
            percent.pop(0)
            sys.stdout.flush()

      t1 = time.time()
      time1 = velocity_verlet.update(dt, Ft, L, eps, sigma, rcutoff, x, p, useCutoff)
      t2 = time.time()
         
      if usePBC:
         boundaries.pbc(L, x, p)     # do the PBC
            
      # scale velocities
      Temp = velocity_verlet.compute_temp(N, m, p)
      LJ.energy(x, p, m, eps, sigma, KE, PE)
      pos[i] = x
      vel[i][:,0] = p[:,0]*scaleP
      vel[i][:,1] = p[:,1]*scaleP
      vel[i][:,2] = p[:,2]*scaleP
      totalKE[i] = np.sum(KE*scaleE)
      totalPE[i] = np.sum(PE*scaleE)
      T[i] = Temp
      t5 = time.time()

      if fp is not None:
         # write current step to file
         np.savetxt(fp, np.c_[pos[i], vel[i], m, KE*scaleE, PE*scaleE, (KE+PE)*scaleE], 
                    fmt="%f %f %f %f %f %f %f %f %f %f", header=str(N)+"\n"+ "x y z vx vy vz m KE PE H", comments = "")

      timeF += t2-t1
      timeE += t5-t2
      
   print "100%. Done!"

   if fp is not None:
#       if skipSteps > 1:
#          # add last step to the file
#          np.savetxt(fp, np.c_[pos[i+1] , vel[i+1], m, totalKE[i+1], totalPE[i+1], totalKE[i+1]+totalPE[i+1]], 
#                     fmt="%f %f %f %f %f %f %f %f %f %f", header=str(N)+"\n"+ "x y z vx vy vz m KE PE H", comments = "")

      fp.close()

   if outFileTimeEvol is not None:
      # save values into file ending with "-eq"
      fpTime = open(outFileTimeEvol, "w")
      fpTime.write("# {0} {1}".format(maxTime, dt))
      np.savetxt(fpTime, np.c_[np.linspace(0,maxTime,M), totalPE, totalKE, totalH, T], fmt="%f %f %f %f %f", header="t PE KE H T")
      fpTime.close()
      print " - energies/temperatures written into: ", outFileTimeEvol

   return totalKE, totalPE, timeF, timeE


def single_run(runCase=1, usePBC=True, useCutoff=False):
   '''
   Main function.
   '''
   runEQ = False         # run equilibrium steps (True) of not (False)
   
   # for Argon; from slogwiki 
   eps = 125.7 * co.Boltzmann
   sigma = 0.3345e-9
   eps = 1
   sigma = 1
   m0 = 1
   #tau = np.sqrt(m*sigma**2/eps)
   Ttarget = 0.002
   Ncalib = 20
   rcutoff = 4.5
   
   if runCase == 3:
      # case 3: specify density
      T = 1
      dt = 1e-3
      N = 100
      n_star = 0.1    # density 
      L = (N*1.0/n_star)**(1/3.)
      m = m0*np.ones(N)
      dist = Distribution()
      x = dist.uniform(N, L, 1.5)
#       p = dist.normal(N, 1e-6)
#       x = dist.lattice(N, L, margin=True)
      p = dist.setvalue(N, 0., 0., 0.)

   elif runCase == 1:
      # case 1: uniform distribution of positions
      T = 1
      dt = 1e-3
      N = 100
      L = 5*N**(1/3.)   # keep the size of the system in relation to number of particles
      m = m0*np.ones(N)
      dist = Distribution()
      x = dist.uniform(N, L, 1.5)
#       p = dist.normal(N, 1e-6)
#       x = dist.lattice(N, L, margin=True)
      p = dist.setvalue(N, 0., 0., 0.)
   
   
   elif runCase == 2:
      # case 2; lattice with minimum distance
      #
      # seems to work with energy increasing after 15 secs. - not conserving the energy!
      #
      T = 20
      dt = 1e-3 # 1e-2*np.sqrt(m0/eps)*sigma/2.
      Nx = 4
      N = Nx**3
      L = 2.21*(Nx-1)
      m = m0*np.ones(N)
      dist = Distribution()
#       p = dist.normal(N, 1e-3)
#       p = dist.setvalue_random_direction(N, 1e-2)
      x = dist.lattice(N, L, margin=0.2)
      p = dist.setvalue(N, 0., 0., 0.)

   else: 
      # toy example: 3 particles
      T = 1e1
      dt = 1e-2
      N = 3
      L = 1.21*(N-1)
      x = np.array([[0.5*L, 0, 0.5*L], [0.5*L, 0.5*L, 0.5*L], [0.5*L, L, 0.5*L]])
      p = np.zeros_like(x)
#       m = np.array([40*co.atomic_mass, 40*co.atomic_mass, 40*co.atomic_mass])
      m = m0*np.ones(N)
#       m[1] *= 2.
      print x, p
   
   print "Number of particles:", N, 
   
   if usePBC:
      print ", using PBC",
   if useCutoff:
      print ", using cutoff radius", rcutoff,
      
   print
   
   maxSteps = 10000
   skipSteps = max(1,int(len(x)/maxSteps))

   outFile = "0_pos_vel_energies.xyz"
   print "Output file (positions, velocities):", outFile, ", data saved every:", skipSteps, " time step(s)"
   outFileEnergies = "0_energies.txt"
   print "Energies/temperature written into: ", outFileEnergies

   params = {'L': L, 'eps': eps, 'sigma': sigma, 'T': T, 'dt': dt, 'Ttarget': Ttarget, 'Ncalibation': Ncalib, 'r_cutoff': rcutoff}

   print params
    
#    plt = Plot()
#    plt.pos_vel_3D(x, p, block=True)
   print "Running case for dt={0}, T={1}, M={2}.".format(dt, T, int(round(T*1.0/dt)))

   t1 = time.time()
   totalKE, totalPE, timeF, timeE = runMD(x, p, m, params, usePBC=usePBC, runEQ=runEQ, useCutoff=useCutoff, outFile=outFile, outFileTimeEvol=outFileEnergies)
   t2 = time.time()
    
#    xLin, pLin, KELin, PELin, HLin = run_case(T, dt, x, p, m, params, LinSpring())
#    xCub, pCub, KECub, PECub, HCub = run_case(T, dt, x, v, m, params, CubicSpring())
#    print "# KE PE H")
#  
   # compute total energies
   totalH = totalKE+totalPE

#    for i in range(10):
#       print i, KE[i], PE[i], H[i])
#       print " 1:", x[i][0], p[i][0])
#       print " 2:", x[i][1], p[i][1])

   print "Energies (min, max): "
   print " - potential: {0}, {1}".format(np.min(totalPE),np.max(totalPE))
   print " - kinetic: {0}, {1}".format(np.min(totalKE),np.max(totalKE))
   print " - total: {0}, {1}".format(np.min(totalH),np.max(totalH))

   timeTotal = t2-t1
   print "Timing:"
   print " - total F: {0:.3f} ({1:.1f} %) ".format(timeF, 100.*timeF/timeTotal)
   print " - total E: {0:.3f} ({1:.1f} %) ".format(timeE, 100.*timeE/timeTotal)
   print " - total time for case: {0:.3f} s".format(timeTotal)
   print
   print "Speed efficiency (time/(M) ): ", (timeTotal/(T*1.0/dt)), "s/time step"


   # write also energies
   outFile = "energies.txt"
   print "Energies written into: ", outFile

   fp = open(outFile, "a")
   fp.write("# {0} {1}".format(T, dt))
   np.savetxt(outFile, np.c_[totalPE, totalKE, totalH], fmt="%f %f %f", header="PE KE H")
   fp.close()

   plt = Plot()
   plt.energies(x, p, T, dt, totalKE, totalPE, totalH, block=False, diff=False);#, fig_base="energies.png")
   plt.energies(x, p, T, dt, totalKE, totalPE, totalH, block=True, diff=True);#, fig_base="diff_energies.png")
#    plt.positions(x, p, T, dt, KE, PE, H)
#    plt.phase_space(x, p, [0,10,20,30,40,50,60,70,80,99], dt, dim='vx-x')
#    plt.pos_vel_3D(x[:,:,0], p[:,:,0], block=False)

#    plt.pos_vel_3D(x[:,:,0:10], p[:,:,0:10])



def timing_runs(runCase=1, num_particles=None, usePBC=True, useCutoff=False):
   '''
   Function to to the timing runs.
   '''
   if num_particles is None:
      num_particles = [10, 10, 100, 1000]
      
   # for Argon; from slogwiki 
   eps = 125.7 * co.Boltzmann
   sigma = 0.3345e-9
   eps = 1
   sigma = 1
   m0 = 1
   tau = 0.5*np.sqrt(m0*sigma**2/eps)
   times = np.zeros((len(num_particles), 3))
   Ttarget = 100
   Ncalib = 100
   rcutoff = 4.5
      
   for iRun, N in enumerate(num_particles):
      # run cases for given number of particles
      if runCase == 1:
         # case 1: uniform distribution of positions
         T = 1
         dt = 1e-3
         L = 5*N**(1/3.)   # keep the size of the system in relation to number of particles
         m = m0*np.ones(N)
         dist = Distribution()
         x = dist.uniform(N, L, 1.5)
#          p = dist.normal(N, 1e-6)
#          x = dist.lattice(N, L, margin=True)
         p = dist.setvalue(N, 0., 0., 0.)
      
      elif runCase == 2:
         # case 2; lattice with minimum distance
         #
         # seems to work with energy increasing after 15 secs. - not conserving the energy!
         #
         T = 1
         dt = 1e-3
         Nx = 4
         N = Nx**3
         L = 2.21*(Nx-1)
         m = m0*np.ones(N)
         dist = Distribution()
#          p = dist.normal(N, 1e-3)
#          p = dist.setvalue_random_direction(N, 1e-2)
         x = dist.lattice(N, L, margin=0.2)
         p = dist.setvalue(N, 0., 0., 0.)
   
      print "Running for # of particles:", N, 
      
      if usePBC:
         print ", using PBC",
      if useCutoff:
         print ", using cutoff radius", rcutoff,
      print
      
#       maxSteps = 10000
#       skipSteps = max(1,int(len(x)/maxSteps))
#       outFile = "0_pos_vel_energies.xyz"
#       print "Output file (positions, velocities):", outFile, ", data saved every:", skipSteps, " time step(s)"
#       outFileEnergies = "energies.txt"
#       print "Energies written into: ", outFileEnergies
   
      params = {'L': L, 'eps': eps, 'sigma': sigma, 'T': T, 'dt': dt, 'Ttarget': Ttarget, 'Ncalibation': Ncalib, 'r_cutoff': rcutoff}
#       print params
    
      print " - running case for dt={0}, T={1}, M={2}.".format(dt, T, int(round(T*1.0/dt)))
      t1 = time.time()
      totalKE, totalPE, timeF, timeE = runMD(x, p, m, params, usePBC=usePBC, runEQ=False, useCutoff=useCutoff)
      t2 = time.time()
      timeTotal = t2-t1
      times[iRun] = [timeF, timeE, timeTotal]

      print " - total time:", timeTotal

   print "Timings done!"

   print times

   
   
   
if __name__ == "__main__":
   # run main function
   if len(sys.argv) == 2 and sys.argv[1].strip().lower() == 'timing':
      # run timing 
      timing_runs(1, useCutoff=True)
   else:
      # "normal" single run
      single_run(1)
   
