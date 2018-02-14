'''
Main module for CMSE 890 group project.

(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek

'''
import time
import sys
import os
import numpy as np
import argparse

import boundaries, input, initialize, integrator, measurables, output

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
   # initialize the system
   pistonPos = params['piston']['z0']        # piston initial position
   pistonVel = params['piston']['v0']        # piston initial velocity
   m = np.full(params['m'])                  # masses
   
   #
   #!! pos, vel = initialize.init(dim, pistonPos, params['n'])
   #
   
   M = int(round(endTime*1.0/dt))         # number of time steps to run

   # time lists for positions, velocities, energies
   posHist = np.zeros((M,N,3))
   velHist = np.zeros_like(posHist)
   KEhist = np.zeros((M, N))
   PEhist = np.zeros_like(KE)
   Ehist = np.zeros(M)
   # progress indicator
   progressList = [int(el) for el in M*np.linspace(0,1.,21)]
   progressList.pop(0)

   if progress:
      print
      print "==> Simulation run"
      print " - computing...",
      
   for i in xrange(1, M):
      # main loop that goes over till given end time
      t = dt*i
      
      if progress:
         # display progress
         if( i == progressList[0] ):
            print "{0:d}%...".format(int((i)*100.0/M)),
            progressList.pop(0)
            sys.stdout.flush()

      # lets do some work here...
      # pos, vel = integrator.velocity_verlet(pos, vel, dt, dim, eps, sigma)

      # impose the boundaries
      pos = boundaries.periodic_boundary_position(pos, N, Lx, Ly)
      pistonPos = boundaries.Piston_Position(pistonPos, pistonVel, dt)
      pos, vel = boundaries.Momentum_Mirror(pos, vel, pistonVel, pistonPos, Lz, N)

      # compute measurables
      #P = measurables.pressure(N, pos, vel, m, dim)
      #T = measurables.temperature(N, pos, vel, m, dim)

      # save values in time history lists
      posHist[i] = pos
      velHist[i] = vel
      KEhist[i] = KE
      PEhist[i] = PE
      Ehist[i] = E

      # output values (?) 
      #??? into which file
      #output.write_step(N, pos, vel, KE, PE, E, T, P)
        
      # prepare for next time step
      # === fill in if necessary ===
  
  
   print "100%. Done!"

   # output for visualization
   #output.write_4_movie()
   #??? we have to define an output file
   #output.write_all_end(posHist, velHist, KEhist, PEhist, Ehist)

   return



if __name__ == '__main__':
   # called from CLI; main function
   args = parse_args()
   
   if args.file:
      # read input file and return read parameters
      params = input.readfile(args.file, args);
   else:
      print("Input file is missing.")
      sys.exit()
      
   runSimulation(params)
