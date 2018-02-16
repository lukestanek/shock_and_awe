'''
Module for parsing input file and checking input parameters. Also defining initial values of the system: dimensions, velocities, ....

(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek
'''
import numpy as np
from numba import jit


def readfile(file, args):
   '''
   Read input file given by filename and other arguments. If same argument if given from command line, use that one, otherwise value from input file is returned.
   
   Args:
      file (string): input file
      args (object): other arguments/parameters given from command line
      
   Returns:
      (dict): system parameters -  
                     params = { 'N': [1, 1, 1],                         # number of particles
                                'spacing': [1.0, 1.0, 1.0],            # initial spacing of particles
                                'm': 0.0,                               # mass of particle
                                'dt': 1e-3, 'endTime': 10,              # time step, end time
                                'piston': {'z0': 0.0, 'v0': 1.0},       # piston init data
                                'eps': 1.0, 'sigma': 1.0,               # LennardJones potential parameters
                                'radius': 0.1                           # cutoff radius used in force calculation
                              }
   '''
   # define initialized dictionary or system parameters
   params = { 'N': [1, 1, 1],                         # number of particles
              'spacing': [1.0, 1.0, 1.0],            # initial spacing of particles
              'm': 0.0,                               # mass of particle
              'dt': 1e-3, 'endTime': 10,              # time step, end time
              'piston': {'z0': 0.0, 'v0': 1.0},       # piston init data
              'eps': 1.0, 'sigma': 1.0,               # LennardJones potential parameters
              'radius': 0.1                           # cutoff radius used in force calculation
            }
   # read file into list ot lines
   fh = open(file, "r")
   lines = fh.readlines()
   fh.close()
   
   for i, line in enumerate(lines):
      # parse all lines
      print(line)
      if line[0] is '#':
         # skip comments
         pass
      
      elif "=" in line:
         #  extract values from file
         fld, restLine = [el.strip() for el in line.split("=")]

         if fld == 'piston':
            # piston has 2 values: z0, v0
            vals = [float(el.strip()) for el in restLine.split(",")]
            params['piston']['z0'] = vals[0]
            params['piston']['v0'] = vals[1]
            
         elif fld == 'N' or fld == 'spacing':
            # number of particles in the system or spacing have 3 values: Nx, Ny, Nz / dx, dy, dz
            vals = [float(el.strip()) for el in restLine.split(",")]
            params[fld] = [ float(vals[0]), float(vals[1]), float(vals[2]) ]

         elif fld in ['m', 'dt', 'endTime', 'eps', 'sigma', 'radius']:
            # parameters that take only one value
            params[fld] = float(restLine.strip())
         
            if type(args) == object and args.getattr(fld) != -1:
               # use value given in CLI
               params[fld] = float(args.getattr(fld))
         else:
            print("Unknown parameter {0} in line {1}".format(line, i))
         
      else: 
         # skip it
         pass

   return params


