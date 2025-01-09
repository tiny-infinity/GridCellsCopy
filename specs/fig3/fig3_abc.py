import numpy as np
"""
LONG SIMULATION TIME~ ~ 6 mins in 40 cores.
This simulation is split into 1000ms segments due to large memory requirements.
"""
def generate_input_params():
    sim_dur = 30000
    input_params = {        
        "sim_dur": sim_dur,
        "sim_id": "VALD-MPD-RAMP-S-s-1a",
        "split_sim":[True,1000],
        "intrnrn_init_noise":[100,0,0.5],
        "stell_init_noise":[100,0,0.5],
        "stell_noise":[sim_dur,0,0],
        "intrnrn_noise":[sim_dur,0,0],
        "stell_const_dc":[4e-4,-2.75e-3],
        "recorder_dt":0.025,
        "init_noise_seed":np.random.randint(100000),
        "record_handle_stell":{"stell_syn_inhib_g": {"state": True,"cells_to_record":"all"},
                               "stell_v":{"state": True,"cells_to_record":"all"}},
        "n_cpus":40
    }
    return input_params
