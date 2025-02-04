"""
Template specs file to run 1D simulation.
"""
import numpy as np

def generate_input_params()-> dict:
    """Generates a dictionary of input parameters for a simulation.
    
    """
    sim_dur=2000 #in ms
    input_params = {   
        "sim_dur": sim_dur,
        "sim_id": "BaseModel", #Important. This is used to save the simulation data
        "init_noise_seed":np.random.randint(0,100000),
        "stell_const_dc":[-2e-3,-7e-3],
        "show_progress_bar":True,
        "n_cpus":8, #all avail
        "split_sim":[False,1000]

    }
    return input_params
