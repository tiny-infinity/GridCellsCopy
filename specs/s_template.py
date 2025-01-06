import numpy as np

def generate_input_params()-> dict:
    """
    Generates a dictionary of input parameters for a simulation.
    Returns:
        dict: A dictionary containing the provided input parameters:
    """
    sim_dur=1000 #in ms
    input_params = {
        "sim_dur": 500,
        "N_intrnrn": 196,
        "N_stell": int(196*2),
        "N_per_sheet": 196, #or ring in 1D
        "sim_id": "BaseModel", #Important. This is used to save the simulation data
        "vel_type": "const",
        "progress":True, #show progress bar. (only for single simulations)
        "init_noise_seed":1000,
        "noise_seed":500,
        "save_conn_matrix":True,
        "vel_type": "input",
        "traj_id": "502",
        "record_handle_stell": {"stell_v": {"state": True,"cells_to_record": [0]}},

    }
    return input_params