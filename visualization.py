'''
Visualization module for CMSE 890 group project.
(c) 2018 - Tom Dixon, Janez Krek, Devin Lake, Ryan Marcus, Luke Stanek
'''

# Import Modules
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

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
    plt.show()
    
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

