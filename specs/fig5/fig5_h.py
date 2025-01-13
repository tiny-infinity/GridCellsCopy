import numpy as np

def generate_input_params():
    sim_dur = float(80000)
    asym_fact=0.5
    is_asym_fact=[1-asym_fact,1+asym_fact] #predictive
    input_params = {
        "sim_dur": sim_dur,
        "sim_id": "VALD-PRED-NET-S-s-1a",
        "vel_type": "PRED-IHD",
        "save_conn_matrix":True,
        "matrix_id":"prosp",
        "g_h_bar":0,
        "intrnrn_init_noise":[100,0,0.5],
        "stell_init_noise":[100,0,0.5],
        "stell_const_dc":[1.4e-3,-2.75e-3],
        "stell_noise": [sim_dur,0,0],
        "intrnrn_noise": [sim_dur,0,0],
        "allothetic_dur":3000,
        "conn_id":"asym", 
        "init_noise_seed":np.random.randint(20000),
        "si_asym_factor":[1,1],
        "is_asym_factor":is_asym_fact,
        "is_peak_asym_fact":1,
        "is_mean_asym_factor":0,
        "recorder_dt":0.25,
        "extra_params":{"stell_dc":0.0019,"dir_change_t":40000},
        "record_handle_stell":{"stell_syn_inhib_g":{"state": True,"cells_to_record":"all"}}

    }
    return input_params