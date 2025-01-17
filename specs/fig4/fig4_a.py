import numpy as np

def generate_input_params():
    sim_dur = float(12000)
    input_params = {
        "sim_dur": sim_dur,
        "show_progress_bar":True,
        "sim_id": "VALD-HCN-SHRK-S-s-3a",
        "intrnrn_init_noise":[100,0,0.5],
        "stell_init_noise":[100,0,0.5],
        "intrnrn_noise":[sim_dur,0,2e-3],
        "stell_noise":[sim_dur,0,1e-3],
        "stell_const_dc":[0.0015,-2.75e-3],
        "init_noise_seed":np.random.randint(0,100000),
        "noise_seed":np.random.randint(0,100000),
        "record_handle_intrnrn":{"intrnrn_v": {"state": False,"cells_to_record":"all"}},
        "record_handle_stell":{"stell_v": {"state": False,"cells_to_record":"all"}}
    }
    return input_params

    #54117-4
    #4125-5
    #1137-1