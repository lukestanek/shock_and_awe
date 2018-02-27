#Import Modules
import numpy as np

# Veloctiy scaling for thermostat
def thermostat(mom, current_temp, desired_temp, N_scale):
    '''
		Objective: Scale the momenta of based on the current temperature 

		Input: Momenta (Nx3 array), current_temp (scalar), desired_temp (scalar), N_scale (scalar)

		Output: scaled momenta (Nx3 array)
	'''

    #print("The temp before scaling is:", np.sum(0.5*px**2 + 0.5*py**2 + 0.5*pz**2))
    
    # Determine momenta scaling
    Scale = np.sqrt(1 + 1/N_scale*(desired_temp/current_temp - 1))
    # Rescale momenta
    mom[:,0] = Scale*mom[:,0]
    mom[:,1] = Scale*mom[:,1]
    mom[:,2] = Scale*mom[:,2]
    
    #print("The temp after scaling is:", np.sum(0.5*px**2 + 0.5*py**2 + 0.5*pz**2))

    return mom
