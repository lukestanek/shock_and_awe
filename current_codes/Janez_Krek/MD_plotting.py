'''

'''
import numpy as np
import sys
import copy


class Plot(object): 
   def results(self, x, p, T, dt, KE, PE, H, fig_base = None, block=True, format='pdf') :
      '''
      Plot results in given arrays and for given end time T and dt.
      
      Args:
         x (dict): position data for all 4 methods - returned from run_case()
         p (dict): momentum data for all 4 methods - returned from run_case()
         T (double): end time
         dt (double): time step
         KE (dict): kinetic energy
         PE (dict): potential energy
         H (dict): total energy
         fig_base (string): format for saving images, {0} is replaced by plot type (angle, velocity...); can include directory (which have to exist before saving file)
      '''
      xAxis = np.arange(0, T, dt)
      N = len(x['fe'][0])
      M = len(xAxis)
   
      import matplotlib.pyplot as plt
   
   #    plt.figure()
   #    plt.grid(True) 
   #    plt.title("forwardEuler: Positions vs. time ($\Delta t=${0}s)".format(dt))
   #    plt.plot(xAxis, [x['fe'][i][0] for i in range(M)], linestyle='dashed', label="wall")
   # 
   #    for j in range(1, N-1):
   #       plt.plot(xAxis, [x['fe'][i][j] for i in range(M)], linestyle='solid', label="mass {0}".format(j))
   #    
   #    plt.plot(xAxis, [x['fe'][i][-1] for i in range(M)], linestyle='dashed', label="wall")
   #    plt.xlabel("time [s]")
   #    plt.ylabel("positions [m]")
   #    plt.legend(loc='upper right')
   #    plt.tight_layout()
   #    
   #    if fig_base != None:
   #       plt.savefig(fig_base.format("positions_fe"), format="png")
   #    else:
   #       plt.show(block=False)
   # 
   #    plt.figure()
   #    plt.grid(True) 
   #    plt.title("Euler-Cromer: Positions vs. time ($\Delta t=${0}s)".format(dt))
   #    plt.plot(xAxis, [x['ec'][i][0] for i in range(M)], linestyle='dashed', label="wall")
   # 
   #    for j in range(1, N-1):
   #       plt.plot(xAxis, [x['ec'][i][j] for i in range(M)], linestyle='solid', label="mass {0}".format(j))
   #    
   #    plt.plot(xAxis, [x['fe'][i][-1] for i in range(M)], linestyle='dashed', label="wall")
   #    plt.xlabel("time [s]")
   #    plt.ylabel("positions [m]")
   #    plt.legend(loc='upper right')
   #    plt.tight_layout()
   #    
   #    if fig_base != None:
   #       plt.savefig(fig_base.format("positions_ec"), format="png")
   #    else:
   #       plt.show(block=False)
   
   
   
   #       plt.figure()
   #       plt.grid(True) 
   #       plt.title("Positions vs. time ($\Delta t=${0:g}s)".format(dt))
   #       plt.plot(xAxis, x['fe'], linestyle='solid', label="forward Euler")
   #       plt.plot(xAxis, x['be'], linestyle='solid', label="backward Euler")
   #       plt.plot(xAxis, x['ec'], linestyle='solid', label="Euler-Cromer")
   #       plt.plot(xAxis, x['an'], linestyle='dashed', label="analytic solution")
   # #       maxX = -1.0 + 1.5*np.max( [x['be'], x['ec'], x['an']] )
   # #       plt.ylim(1.0-maxX, 1.0+maxX)
   #       plt.xlabel("time [s]")
   #       plt.ylabel("position [m]")
   #       plt.legend(loc='lower left')
   #       plt.tight_layout()
        
   #       if fig_base != None:
   #          plt.savefig(fig_base.format("angle"), format="png")
   #       else:
   #          plt.show(block=False)
   
   #       plt.show()
   
      return 0
   
   
   def positions(self, x, p, T, dt, KE, PE, H, dim='x', fig_base = None, block=True, format='pdf') :
      '''
      Plot results in given arrays and for given end time T and dt.
      
      Args:
         x (dict): position data for all 4 methods - returned from run_case()
         p (dict): momentum data for all 4 methods - returned from run_case()
         T (double): end time
         dt (double): time step
         KE (dict): kinetic energy
         PE (dict): potential energy
         H (dict): total energy
         dim (string): dimension to plot: x, y, z
         fig_base (string): format for saving images, {0} is replaced by plot type (angle, velocity...); can include directory (which have to exist before saving file)
      '''
      xAxis = np.arange(0, T, dt)
      N = len(x[0])
      M = len(xAxis)
   
      import matplotlib.pyplot as plt
   
      plt.figure()
      plt.grid(True) 
      plt.title("Positions vs. time ($\Delta t=${0}s)".format(dt))
    
      for j in range(1, N-1):
         plt.plot(xAxis, x[:,j], linestyle='solid', label="part {0}".format(j))
       
      plt.xlabel("time [s]")
      plt.ylabel("positions [m]")
      plt.legend(loc='upper right')
      plt.tight_layout()
       
      if fig_base != None:
         plt.savefig(fig_base.format("positions"), format=format)
      else:
         plt.show(block)
    
      return 0
   
   
   def energies(self, x, p, T, dt, KE, PE, H, fig_base = None, block=True, format='pdf', diff=False) :
      '''
      Plot results in given arrays and for given end time T and dt.
      
      Args:
         x (dict): position data for all 4 methods - returned from run_case()
         p (dict): momentum data for all 4 methods - returned from run_case()
         T (double): end time
         dt (double): time step
         KE (dict): kinetic energy
         PE (dict): potential energy
         H (dict): total energy
         fig_base (string): format for saving images, {0} is replaced by plot type (angle, velocity...); can include directory (which have to exist before saving file)
      '''
      xAxis = np.arange(0, T-0.5*dt, dt)
      N = len(x[0])
      M = len(xAxis)
   
      import matplotlib.pyplot as plt
#       plt.figure(figsize=(20, 6), dpi=300)
      plt.figure()
      plt.grid(True) 
      plt.title("Energy vs. time ($\Delta t=${0}s)".format(dt))
   
      if diff is True:
         # plot difference between total energies
         plt.plot(xAxis, H/H[0] - 1, linestyle='solid', label="diff")
         plt.ylabel("difference in total energy []")
         plt.legend(loc='upper right')
      else:
         # plot all energies in absolute values
         plt.plot(xAxis, KE, linestyle='solid', label="kinetic")
         plt.plot(xAxis, PE, linestyle='solid', label="potential")
         plt.plot(xAxis, H, linestyle='solid', label="total")
         plt.ylabel("energy [J]")
         plt.legend(loc='center right')

      plt.xlabel("time [s]")
      plt.tight_layout()
   
      if fig_base != None:
         plt.savefig(fig_base.format("energies"), format=format)
      else:
         plt.show(block)
   
      return 0
   
   
   def phase_space(self, x, p, ti, dt, dim='vx-x', fig_base = None, block=True, format='pdf') :
      '''
      Plot phase-space as movie in given arrays and for given end time T and dt.
      
      Args:
         x (dict): position data for all 4 methods - returned from run_case()
         p (dict): momentum data for all 4 methods - returned from run_case()
         ti (int or list of integers): time index to plot
         dt (double): time step
         KE (dict): kinetic energy
         PE (dict): potential energy
         H (dict): total energy
         dim (string): which phase space to plot - velocity and coordinate are separated by dash '-'; vx-x, vy-x, ...
         fig_base (string): format for saving images, {0} is replaced by plot type (angle, velocity...); can include directory (which have to exist before saving file)
      '''
      N = len(x[0])
      # check given coordinates and get coordinates and their indexes
      coordList = ['x', 'y', 'z']
      vcoord, xcoord = dim.split("-")
      
      xcoord = "x" if xcoord not in coordList else xcoord
      vcoord = "x" if vcoord not in coordList else vcoord
      xi = coordList.index(xcoord)
      vi = coordList.index(vcoord)
   
      import matplotlib.pyplot as plt
      # plt.figure(figsize=(10, 14), dpi=300)
      fig = plt.figure()
      plt.grid(True) 
      
      if type(ti) is not list:
         plt.title("v{0}-{1} phase-space at t={2} ($\Delta t=${3}s)".format(vcoord, xcoord, ti*dt, dt))
         plt.scatter(x[ti][:,xi], p[ti][:,vi])
         plt.tight_layout()
      else:
         plt.title("v{0}-{1} phase-space ($\Delta t=${2}s)".format(vcoord, xcoord, dt))
         
         for j in ti:
            plt.scatter(x[j][:,xi], p[j][:,vi], label="t={0:g}s".format(j*dt))

         plt.legend(bbox_to_anchor=(1.01,1.0))
         plt.tight_layout(rect=(0.05,0.05,0.85,1))

      plt.xlabel("positions [m]")
      plt.ylabel("momentum [Nm]")

      if fig_base != None:
         plt.savefig(fig_base.format("phase_space_time"), format=format)
      else:
         plt.show(block)
   
      return 0
 
 
   def phase_space_time(self, x, p, T, dt, KE, PE, H, dim='vx-x', fig_base = None, block=True, format='pdf') :
      '''
      Plot phase-space as movie in given arrays and for given end time T and dt.
      
      Args:
         x (dict): position data for all 4 methods - returned from run_case()
         p (dict): momentum data for all 4 methods - returned from run_case()
         T (double): end time
         dt (double): time step
         KE (dict): kinetic energy
         PE (dict): potential energy
         H (dict): total energy
         dim (string): which phase space to plot - velocity and coordinate are separated by dash '-'; vx-x, vy-x, ...
         fig_base (string): format for saving images, {0} is replaced by plot type (angle, velocity...); can include directory (which have to exist before saving file)
      '''
      xAxis = np.arange(0, T, dt)
      N = len(x[0])
      M = len(xAxis)
      pList = [int(round(el)) for el in np.linspace(1, N-2, 10)]
      # check given coordinates and get coordinates and their indexes
      coordList = ['x', 'y', 'z']
      vcoord, xcoord = dim.split("-")
      
      xcoord = "x" if xcoord not in coordList else xcoord
      vcoord = "x" if vcoord not in coordList else vcoord
      xi = coordList.index(xcoord)
      vi = coordList.index(vcoord)
   
      import matplotlib.pyplot as plt
      # plt.figure(figsize=(10, 14), dpi=300)
      plt.figure()
      plt.grid(True) 
      plt.title("{0}-{1} phase-space of particles: {2} ($\Delta t=${3}s)".format(vcoord, xcoord, ", ".join([str(el) for el in pList]), dt))
    
      for pj in range(len(pList)):
         j = pList[pj]
         plt.plot(x[:,j][xi], p[:,j][vi], linestyle='solid', linewidth=1)
       
      plt.xlabel("positions [m]")
      plt.ylabel("momentum [Nm]")
      plt.tight_layout()
       
      if fig_base != None:
         plt.savefig(fig_base.format("phase_space_time"), format=format)
      else:
         plt.show(block)
   
      return 0

   

   def pos_vel_3D(self, x, v, block=True, fig_base=None, format='pdf'):
      import matplotlib.pyplot as plt
      from mpl_toolkits.mplot3d import Axes3D
      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')
      ax.scatter(x[:,0], x[:,1], x[:,2], c='red')
      ax.set_xlabel('x')
      ax.set_ylabel('y')
      ax.set_zlabel('z')
      plt.title("Positions")
      plt.tight_layout()
      plt.show(block=False)
   
      import matplotlib.pyplot as plt
      from mpl_toolkits.mplot3d import Axes3D
      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')
      ax.scatter(v[:,0], v[:,1], v[:,2], c='blue')
      ax.set_xlabel('$v_x$')
      ax.set_ylabel('$v_y$')
      ax.set_zlabel('$v_z$')
      plt.title("Velocities")
      plt.tight_layout()

      if fig_base != None:
         plt.savefig(fig_base.format("phase_space"), format=format)
      else:
         plt.show(block)
     
      plt.show(block=False)

 
 
 
def plot_temperasture(self):
   d = np.loadtxt('0_energies.txt-eq')
   d1 = np.loadtxt('0_energies.txt')
   
   plt.grid(True);
   plt.title('Temperature vs. time')
   plt.plot(d[:,0], d[:,4],label='w/ thermostat')
   plt.plot(d1[:,0], d1[:,4],label='w/o thermostat')
   plt.xlabel('time [s]')
   plt.ylabel('temperature [eV]')
   plt.legend()
   plt.tight_layout()
   plt.show()
    
    
def plot_with_without_cutoff():
   Nlist = [10, 100, 1000]
   Msteps = 1000
   
   timesNoCutoff =  np.array([[6.47962093e-02, 6.61993027e-02, 1.34366035e-01],
                               [7.04482341e+00, 1.33016014e+00, 8.39151692e+00],
                               [7.37698497e+02, 1.25663543e+02, 8.64312157e+02]])
   timesWithCutoff =  np.array([[5.80692291e-02, 1.36533737e-01, 2.01596975e-01],
                               [2.80293536e+00, 1.51875567e+00, 4.33761907e+00],
                               [2.41038569e+02, 1.23165506e+02, 3.64691274e+02]])
   timesNoCutoff /= Msteps
   timesWithCutoff /= Msteps
   
   fig = plt.figure()
   plt.grid(True);
   plt.title('Computing time of forces vs. number of particles')
   plt.loglog(Nlist, timesNoCutoff[:,0],label='w/o cutoff radius')
   plt.loglog(Nlist, timesWithCutoff[:,0],label='w/ cutoff radius 4.5$\sigma$')
   plt.xlabel('# of particles')
   plt.ylabel('time [s]')
   plt.legend()
   plt.tight_layout()
   plt.show(block=False)
    
   fig = plt.figure()
   plt.grid(True);
   plt.title('Total computing time vs. number of particles')
   plt.loglog(Nlist, timesNoCutoff[:,2],label='w/o cutoff radius')
   plt.loglog(Nlist, timesWithCutoff[:,2],label='w/ cutoff radius 4.5$\sigma$')
   plt.xlabel('# of particles')
   plt.ylabel('time [s]')
   plt.legend()
   plt.tight_layout()
   plt.show()
    
    
    
    
    
    