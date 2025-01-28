"""
Template specs file to run multiple simulations.
"""

import numpy as np

def generate_mult_input_params()-> dict:
    """Generates a dictionary of multiple input parameters

    Used for multiple simulations. 
    Each input_params dictionary should contain a unique sim_num.
    
    Returns:
        dict: A dictionary all input parameters
    """
    sim_dur=200 #in ms
    mult_input_params = {}
    for sim_num in range(3):
        input_params = {
            "sim_dur": sim_dur,
            "N_intrnrn": 196,
            "N_stell": int(196*2),
            "N_per_sheet": 196, #or ring in 1D
            
            "sim_id": "mBaseModel",  #Important. 
                                    #Used to save the simulation data
            
            "sim_num":str(sim_num)       #Important.
                                    # Assign a simulation number to 
                                    # every simulation

        }
        mult_input_params[str(sim_num)] = input_params #Important.
    return mult_input_params