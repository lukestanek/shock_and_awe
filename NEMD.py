'''
Main module for CMSE 890 group project.

(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek

'''
import time
import sys
import os
import numpy as np
import argparse

import boundaries, input, initialization, integrator, measurables, output

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
   Lx, Ly, Lz = params['Lx'], params['Lz'], params['Lz']
   dim = np.array([[0,Lx], [0,Ly], [0,Ly]])      # x, y, z
   eps = params['eps']
   sigma = params['sigma']
   radius = params['radius']
   # initialize the system
   pistonPos = params['piston']['z0']        # piston initial position
   pistonVel = params['piston']['v0']        # piston initial velocity
   
   #
   #!! pos, mom = initialize.init(dim, pistonPos, params['n'])
   #
   
   print(params)
   N = 10
   
   pos = np.zeros((N,3))

   # assign positions in all three coordinates
   pos[:,0] = Lx * np.random.random(N)
   pos[:,1] = Ly * np.random.random(N)
   pos[:,2] = Lz * np.random.random(N)


   print(pos)

   mom = np.full((N,3),0)
   
   
   
   
   M = int(round(endTime*1.0/dt))            # number of time steps to run
   m = np.full(N, params['m'])                  # masses

   # time lists for positions, velocities, energies
   posHist = np.zeros((M,N,3))
   momHist = np.zeros_like(posHist)
   KEhist = np.zeros((M, N))
   PEhist = np.zeros_like(KEhist)
   Ehist = np.zeros(M)
   # progress indicator
   progressList = [int(el) for el in M*np.linspace(0,1.,21)]
   progressList.pop(0)

   if progress:
      print("==> Simulation run\n - computing...",)
      
   # compute forces on initial particles
   force = integrator.calc_force(pos, radius, Lx, Ly)
   
   for i in range(1, M):
      # main loop that goes over till given end time
      t = dt*i
      print(t, end='') ; sys.stdout.flush()
      
      if progress:
         # display progress
         if( i == progressList[0] ):
            print("{0:d}%...".format(int((i)*100.0/M)),)
            progressList.pop(0)
            sys.stdout.flush()

      # lets do some work here...
      print(" 1", end='') ; sys.stdout.flush()
      # update positions+momentum
      pos, mom, force = integrator.vel_ver(pos, mom, pistonPos, pistonVel, dt, force, Lx, Ly, Lz, radius)
      # impose the boundaries
      print(" 2", end='') ; sys.stdout.flush()
      pos = boundaries.periodic_boundary_position(pos, N, Lx, Ly)
      print(" 3", end='') ; sys.stdout.flush()
      pistonPos = boundaries.Piston_Position(pistonPos, pistonVel, dt)
      print(" 4", end='') ; sys.stdout.flush()
      pos, mom = boundaries.Momentum_Mirror(pos, mom, pistonVel, pistonPos, Lz, N)
      print(" 5", end='') ; sys.stdout.flush()

      # compute measurables
      #P = measurables.pressure(N, pos, mom, m, dim)
      #T = measurables.temperature(N, pos, mom, m, dim)

      # save values in time history lists
      posHist[i] = pos
      momHist[i] = mom
      KEhist[i] = KE
      PEhist[i] = PE
      Ehist[i] = E

      # output values (?) 
      #??? into which file
      #output.write_step(N, pos, mom, KE, PE, E, T, P)
        
      # prepare for next time step
      # === fill in if necessary ===
      print("")
  
   print("100%. Done!")

   # output for visualization
   #output.write_4_movie()
   #??? we have to define an output file
   #output.write_all_end(posHist, momHist, KEhist, PEhist, Ehist)

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
