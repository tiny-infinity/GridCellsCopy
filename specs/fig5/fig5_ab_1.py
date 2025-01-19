import numpy as np
def generate_input_params():
    sim_dur = float(80000)
    input_params = {
        "sim_dur": sim_dur,
        "sim_id": "VALD-PRED-INT-S-s-1a",
        "vel_type": "PRED-IHD",
        "save_conn_matrix":True,
        "split_sim":[True,8000],
        "matrix_id":"default",
        "g_h_bar":0.0015,
        "intrnrn_init_noise":[100,0,0.5],
        "stell_init_noise":[100,0,0.5],
        "stell_const_dc":[1.4e-3,-2.75e-3],
        "stell_noise": [sim_dur,0,0],
        "intrnrn_noise": [sim_dur,0,0],
        "init_noise_seed":100,
        "n_phases":64,
        "recorder_dt":0.1,
        "hs_tau": 5.6,
        "hf_tau": 0.51, 
        "conn_id":"1d",
        "allothetic_dur":3000, 
        "si_asym_factor":[1,1],
        "is_asym_factor":[1,1],
        "is_peak_asym_fact":1,
        "is_mean_asym_factor":0,
        "extra_params":{"stell_dc":-2e-3,"dir_change_t":40000},
        "record_handle_stell":{"stell_syn_inhib_g":{"state": True,"cells_to_record":"all"}}
    }
    return input_params