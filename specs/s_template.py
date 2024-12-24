import numpy as np

#Template to run a single 1D simulation using s_sim_setup.py
"""
We use a function instead of directly defining a dictionary to allow for precalculations to set cetain parameters. 
For e.g if one parameter is a function of another, 
    we can define the function in the generate_input_params function and use it to set the value of the parameter.
"""
def generate_input_params()-> dict:
    """
    Generates a dictionary of input parameters for a simulation.
    Returns:
        dict: A dictionary containing the provided input parameters:
    """
    sim_dur=1000 #in ms
    input_params = {
        "sim_dur": sim_dur,
        "n_intrnrn": 196,
        "n_stell": int(196*2),
        "N_per_sheet": 196, #or ring in 1D
        "sim_id": "BaseModel", #Important. This is used to save the simulation data
        "vel_type": "const",
        "progress":False, #show progress bar. (only for single simulations)
        "init_noise_seed":1000,
        "noise_seed":500,
        "save_conn_matrix":True,
    }
    return input_params