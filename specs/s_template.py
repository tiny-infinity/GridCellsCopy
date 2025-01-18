import numpy as np

def generate_input_params()-> dict:
    """
    Generates a dictionary of input parameters for a simulation.
    Returns:
        dict: A dictionary containing the provided input parameters:
    """
    sim_dur=20000 #in ms
    input_params = {   
        "sim_dur": sim_dur,
        "sim_id": "BaseModel", #Important. This is used to save the simulation data
        "init_noise_seed":np.random.randint(0,100000),
        "stell_const_dc":[-2e-3,-7e-3],

    }
    return input_params
