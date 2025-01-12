import numpy as np
def generate_input_params():
    sim_dur = float(40000)

    input_params = {
        "sim_dur": sim_dur,
        "sim_id": "VALD-HCN-SHRK-S-s-2a", #"VI-TI-SIMS-S-1a",
        "si_peak":0,
        "g_h_bar":0,
        "intrnrn_init_noise":[100,0,0.5],
        "stell_init_noise":[100,0,0.5],
        "stell_const_dc":[0.0015,-1],
        "init_noise_seed":50,
        "record_handle_stell":{"stell_v":{"state": True,"cells_to_record":"all"}}
    }
    return input_params
    #54117-4