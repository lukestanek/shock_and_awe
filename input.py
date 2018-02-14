'''
Module for parsing input file and checking input parameters. Also defining initial values of the system: dimensions, velocities, ....

(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek
'''

def readfile(file, args):
   '''
   Read input file given by filename and other arguments. If same argument if given from command line, use that one, otherwise value from input file is returned.
   
   Args:
      file (string): input file
      args (object): other arguments/parameters given from command line
      
   Returns:
      (dict): system parameters -  
               params = {'Lx': 1.0, 'Ly': 1.0, 'Lz': 1.0,         # system dimension
                          'n': 1.0,                               # initial particle density
                          'm': 0.0,                               # mass of particle
                          'dt': 1e-3, 'endTime': 10,              # time step, end time
                          'piston': {'z0': 0.0, 'v0': 1.0},       # piston init data
                          'eps': 1.0, 'sigma': 1.0,               # LennardJones potential parameters
                        }
   '''
   # define initialized dictionary or system parameters
   params = {'Lx': 1.0, 'Ly': 1.0, 'Lz': 1.0,         # system dimension
              'n': 1.0,                               # initial particle density
              'm': 0.0,                               # mass of particle
              'dt': 1e-3, 'endTime': 10,              # time step, end time
              'piston': {'z0': 0.0, 'v0': 1.0},       # piston init data
              'eps': 1.0, 'sigma': 1.0,               # LennardJones potential parameters
            }
   # read file into list ot lines
   fh = open(file, "r")
   lines = fh.readlines()
   fh.close()
   
   for i, line in enumerate(lines):
      # parse all lines
      print(line)
      if "=" in line:
         #  extract values from file
         fld, restLine = [el.strip() for el in line.split("=")]

         if fld == 'piston':
            # piston has 2 values: z0, v0
            vals = [float(el.strip()) for el in restLine.split(",")]
            params['piston']['z0'] = vals[0]
            params['piston']['v0'] = vals[1]
            
         elif fld == 'size':
            # dimension of the system has 3 values: Lx, Ly, Lz
            vals = [float(el.strip()) for el in restLine.split(",")]
            params['Lx'] = float(vals[0])
            params['Ly'] = float(vals[1])
            params['Lz'] = float(vals[2])

         elif fld in ['m', 'dt', 'endTime', 'eps', 'sigma', 'n']:
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


