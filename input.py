'''
Module for parsing input file and checking input parameters. Also defining initial values of the system: dimensions, velocities, ....

(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek


History:
v0.1    - JK, 2018-02 initial read of input file
v0.2    - JK, 2018-02-22 -- added end time of piston movement; input filename is part of parameter dict

'''
import numpy as np
from numba import jit
import argparse

def readfile(file, args):
    '''
    Read input file given by filename and other arguments. If same argument if given from command line, use that one, otherwise value from input file is returned.
    
    Args:
       file (string): input file
       args (object): other arguments/parameters given from command line
       
    Returns:
       (dict): system parameters -  
                    params = { 'N': [1, 1, 1],                                     # number of particles
                                'spacing': [1.0, 1.0, 1.0],                         # initial spacing of particles
                                'm': 0.0,                                           # mass of particle
                                'dt': 1e-3, 'endTime': 10,                          # time step, end time
                                'piston': {'z0': 0.0, 'v0': 1.0, 'endT': -1},       # piston init data: init positon, velocity, end time of movement
                                'eps': 1.0, 'sigma': 1.0,                           # LennardJones potential parameters
                                'radius': 0.1,                                      # cutoff radius used in force calculation
                                'input_filename': ''                                # input file name
                                'Tdesired': -1,                                     # desired temperatue in the system; -1 means no thermostat
                                'eq_run': {'M_eq': 0, "M_scale": 0}                 # time steps in equlibration run and ratio for thermostat 
                               }
    '''
    # define initialized dictionary or system parameters
    params = { 'N': [1, 1, 1],                                     # number of particles
               'spacing': [1.0, 1.0, 1.0],                         # initial spacing of particles
               'm': 0.0,                                           # mass of particle
               'dt': 1e-3, 'endTime': 10,                          # time step, end time
               'piston': {'z0': 0.0, 'v0': 1.0, 'endT': -1},       # piston init data: init positon, velocity, end time of movement (-1 = no stopping time)
               'eps': 1.0, 'sigma': 1.0,                           # LennardJones potential parameters
               'radius': 0.1,                                      # cutoff radius used in force calculation
               'input_filename': '',                               # input file name
               'Tdesired': -1,                                     # desired temperatue in the system; -1 means no thermostat
               'eq_run': {'M_eq': 0, "M_scale": 0}                 # time steps in equlibration run and ratio for thermostat 
    }
    # read file into list ot lines
    fh = open(file, "r")
    lines = fh.readlines()
    fh.close()
    
    cli_params = vars(args)
    params['input_filename'] = file
#     import sys
#     sys.exit()
    
    for i, line in enumerate(lines):
        # parse all lines
        if line[0] is '#':
            # skip comments
            pass
            
        elif "=" in line:
            #  extract values from file
            dataInLine = line.split("#")[0]     # remove any comments in the last part of line
            fld, restLine = [el.strip() for el in dataInLine.strip().split("=")]
            print(i, fld)
            
            if fld == 'piston':
                # piston has 3 values: z0, v0, endtime of moving
                vals = [float(el.strip()) for el in restLine.split(",")]
                params['piston']['z0'] = vals[0]
                params['piston']['v0'] = vals[1]
                params['piston']['endT'] = vals[2]
                
            elif fld == 'N':
                #  number of particles in the system has 3 values: Nx, Ny, Nz
                vals = [float(el.strip()) for el in restLine.split(",")]
                params[fld] = [ int(vals[0]), int(vals[1]), int(vals[2]) ]
                
            elif fld == 'spacing':
                # spacing has 3 values: dx, dy, dz
                vals = [float(el.strip()) for el in restLine.split(",")]
                params[fld] = [ float(vals[0]), float(vals[1]), float(vals[2]) ]

            elif fld == 'eq_run':
                # equilibration data: # step for eq. and Nstep for thermostat
                vals = [float(el.strip()) for el in restLine.split(",")]
                params['eq_run']['M_eq'] = int(vals[0])
                params['eq_run']['M_scale'] = float(vals[1])
                
            elif fld in ['m', 'dt', 'endTime', 'eps', 'sigma', 'radius', 'Tdesired']:
                # parameters that take only one value
                params[fld] = float(restLine.strip())
                
            else:
                print("Unknown parameter {0} in line {1}".format(line, i))
        
        else: 
           # skip it
           pass

    if 'eq_run' in cli_params.keys():
        # manage eq_run parames from CLI
        M_eq, M_scale = args.eq_run
        params['eq_run']['M_eq'] = int(M_eq)
        params['eq_run']['M_scale'] = float(M_scale)
        print("Set {0} to M_eq={1}, M_scale={2} (CLI param)".format(fld, params['eq_run']['M_eq'], params['eq_run']['M_scale']))
        
    elif fld in cli_params.keys() and cli_params[fld] != -1:
        # use value given in CLI
        print("Set {0} to {1} (CLI param)".format(fld, cli_params[fld]))
        params[fld] = float(cli_params[fld])

    return params


