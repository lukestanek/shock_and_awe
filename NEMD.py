'''
Main module for CMSE 890 group project.

(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek

History:
v0.1    - JK, 2018-02?? -- initial
v0.2    - JK, 2018-02-22 -- piston end time added;

'''
import time
import sys
import os
import numpy as np
import argparse
import time

import boundaries, input, initialization, integrator, measurables, output, visualization

def parse_args():
   '''
   Define CLI parameters script will accept. Except "file" all parameters are optional.
   
   Returns:
      (object): parsed arguments with given values 
   '''
   argparser = argparse.ArgumentParser()
   argparser.add_argument("-endTime", help="stop simulation at given time", nargs='?', default=-1)
   argparser.add_argument("-dt", help="time step to use in simulation", nargs='?', default=-1)
   argparser.add_argument("file", help="name of input file", nargs=1)
   args = argparser.parse_args()
   return args



def runSimulation(params):
   '''
   Main function that runs the simulation with given parameters.
   
   Args:
      params (dict): parameters for running the simulation (see input.readfile for description of parameters)
      
   '''
   progress = True
   # setup values often used
   t = 0
   endTime = params['endTime']
   dt = params['dt']
   num_particles = params['N']               # number of particles in all three dimensions
   spacing = params['spacing']               # particle spacing in all three dimensions
   eps = params['eps']
   sigma = params['sigma']
   radius = params['radius']
   # initialize the system
   pistonPos = params['piston']['z0']        # piston initial position
   pistonVel = params['piston']['v0']        # piston initial velocity
   pistonEndTime = params['piston']['endT']  # time when piston stops moving

   print(params)
   
   # initialize particle positions and momentums; also get system length
   pos, mom, [Lx, Ly, Lz] = initialization.initilization(num_particles, spacing)
   N = len(pos)
   print("System dimensions (x,y,z): {0}, {1}, {2}".format(Lx, Ly, Lz))
   dim = np.array([[0,Lx], [0,Ly], [0,Lz]])      # x, y, z
   
   M = int(round(endTime*1.0/dt))            # number of time steps to run
   m = np.full(N, params['m'])               # masses

   # time lists for positions, velocities, energies
   pistHist = np.zeros(M)
   posHist = np.zeros((M,N,3))
   momHist = np.zeros_like(posHist)
   KEhist = np.zeros((M, N))
   PEhist = np.zeros_like(KEhist)
   Ehist = np.zeros(M)
   # progress indicator
   progressList = [int(el) for el in M*np.linspace(0,1.,21)]
   progressList.pop(0)

   posHist[0] = pos
   momHist[0] = mom
   tStart = time.time()

   if progress:
      print("==> Simulation run\n - computing...", end='')
      
   # compute forces on initial particles
   force = integrator.calc_force(pos, radius, Lx, Ly)

   for i in range(1, M):
      # main loop that goes over till given end time
      t = dt*i
      
      if progress:
         # display progress
         if( i == progressList[0] ):
            print("{0:d}%...".format(int((i)*100.0/M)),end='', flush=True)
            progressList.pop(0)

      # lets do some work here...
      # update positions+momentum
      pos, mom, force = integrator.vel_ver(pos, mom, pistonPos, pistonVel, dt, force, Lx, Ly, Lz, radius)
      # impose the boundaries
      pos = boundaries.periodic_boundary_position(pos, N, Lx, Ly)
      pistonPos = boundaries.calc_Piston_Position(pistonPos, pistonVel, dt, pistonEndTime, t)
      pos, mom = boundaries.Momentum_Mirror(pos, mom, pistonVel, pistonPos, Lz, dt, N)

      # compute measurables
      KE = measurables.calc_kinetic(mom)
      #P = measurables.pressure(N, pos, mom, m, dim)
      T = measurables.calc_temp(mom)

      # save values in time history lists
      pistHist[i] = pistonPos #z component of pistion
      posHist[i] = pos
      momHist[i] = mom
      KEhist[i] = KE
#       PEhist[i] = PE
#       Ehist[i] = KE+PE

      # output values (?) 
      #??? into which file
      #output.write_step(N, pos, mom, KE, PE, E, T, P)
        
      # prepare for next time step
      # === fill in if necessary ===
  
   print("100%. Done!")

   tEnd = time.time()
   
   print(" - dt = {0:.2f}s ".format(tEnd-tStart))

   # output for visualization
   #output.write_4_movie()
   #??? we have to define an output file
   #output.write_all_end(posHist, momHist, KEhist, PEhist, Ehist)

   # get input file base name to use it as part of output files...
   baseName, fileext = os.path.splitext(params['input_filename'])
   outFile = "0_{0}_pos_vel_KE.txt".format(baseName)
   print(" - writing into file: {0}".format(outFile))
   output.write_pos_vel_hist(outFile, posHist, momHist, KEhist, pistHist, Lx, Ly, Lz)

   # visualize initial and end positions 
   #visualization.visualize(posHist[0],momHist[0])
   #visualization.visualize(posHist[-1],momHist[-1])

   return



if __name__ == '__main__':
   # called from CLI; main function
   args = parse_args()
   
   if args.file:
      # read input file and return read parameters
      params = input.readfile(args.file[0], args);
   else:
      print("Input file is missing.")
      sys.exit()
      
   runSimulation(params)
