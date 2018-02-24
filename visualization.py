'''
Visualization module for CMSE 890 group project.
(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek
'''

# Import Modules
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def visualize(pos,mom):
    
    '''
        Objective: Visualize the current stat of the simulation
        
        Inputs: 
                pos: Position array of the system
                mom: Momenta array of the system
                
        Outpus:
                none
    '''

    # Plotting Position
    fig1=plt.figure(figsize=(8,8))
    ax = plt.axes(projection='3d')
    ax.scatter(pos[:,2], pos[:,0], pos[:,1], c=pos[:,2], cmap='viridis', linewidth=6);
    plt.title('Initialized Position Lattice Like with Perturbation')
    ax.set_xlabel('z')
    ax.set_ylabel('x')
    ax.set_zlabel('y')
    plt.show(block=False)
    
    # Plotting Momenta
    fig2 = plt.figure(figsize=(8,8))
    ax2 = plt.axes(projection='3d')
    ax2.scatter(mom[:,2], mom[:,0], mom[:,1], c=mom[:,2], cmap='viridis', linewidth=6);
    plt.title('Initialized Momenta (Uniform Mean 0 and width 1)')
    ax.set_xlabel('pz')
    ax.set_ylabel('px')
    ax.set_zlabel('py')
    plt.show()
    
    
    return 


def energies(T, dt, KE, PE, H, fig_base = None, block=True, format='pdf', diff=False) :
    '''
    Plot results in given arrays and for given end time T and dt.
    
    Args:
        T (double): end time
        dt (double): time step
        KE (dict): kinetic energy
        PE (dict): potential energy
        H (dict): total energy
        fig_base (string): format for saving images, {0} is replaced by plot type (angle, velocity...); can include directory (which have to exist before saving file)
    '''
    xAxis = np.arange(0, T-0.5*dt, dt)
    M = len(xAxis)
    
    import matplotlib.pyplot as plt
    #       plt.figure(figsize=(20, 6), dpi=300)
    plt.figure()
    plt.grid(True) 
    
    if diff:
        plt.title("Energy difference vs. time ($\Delta t=${0}s)".format(dt))
    else:
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
