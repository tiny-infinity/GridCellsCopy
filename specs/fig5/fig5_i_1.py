import numpy as np

def generate_mult_input_params():
    mult_input_params = {}
    sim_dur = float(80000)
    asym_factor_arr = np.linspace(0.1,0.9,9)
    n_trials = 1
    sim_num = 0
    for asym_fact in asym_factor_arr:
        for tr in range(n_trials):
            input_params = {
                "sim_num":sim_num,
                "sim_dur": sim_dur,
                "sim_id": "VALD-PRED-NET-S-m-1a",
                "vel_type": "PRED-IHD",
                "g_h_bar":0,
                "intrnrn_init_noise":[100,0,0.5],
                "stell_init_noise":[100,0,0.5],
                "stell_const_dc":[1.4e-3,-2.75e-3],
                "stell_noise": [sim_dur,0,0],
                "intrnrn_noise": [sim_dur,0,0],
                "n_phases":64,
                "recorder_dt":0.25,
                "allothetic_dur":3000,
                "hs_tau": 5.6,
                "hf_tau": 0.51,
                "conn_id":"asym",
                "n_nodes":3, 
                "data_root":"data/",
                "init_noise_seed":np.random.randint(100000),
                "si_asym_factor":[1,1],
                "is_asym_factor":[1-asym_fact,1+asym_fact],
                "is_peak_asym_fact":1,
                "is_mean_asym_factor":0,
                "extra_params":{"stell_dc":0.0019,"dir_change_t":40000},
                }
            mult_input_params[sim_num] = input_params
            sim_num +=1
    return mult_input_params
