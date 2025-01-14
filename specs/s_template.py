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
        "vel_type": "const",
        "init_noise_seed":np.random.randint(0,100000),
        "intrnrn_init_noise":[100,0,0.5],
        "stell_init_noise":[100,0,0.5],
        "stell_const_dc":[-7.5e-4,-3e-2],

    }
    return input_params
