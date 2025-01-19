import numpy as np

def generate_mult_input_params():
    sim_dur = float(80000)
    gh_range_arr = np.linspace(0,0.002,21)
    n_trials = 1
    sim_num = 0
    mult_input_params = {}
    for gh in gh_range_arr:
        for tr in range(n_trials):
            input_params = {
                "sim_num":str(sim_num),
                "sim_dur": sim_dur,
                "sim_id": "VALD-PRED-INT-S-m-1b",
                "g_h_bar":gh,
                "n_nodes":3,
                "data_root":"data/",
                "intrnrn_init_noise":[100,0,0.5],
                "stell_init_noise":[100,0,0.5],
                "stell_const_dc":[0.002,-2.75e-3],
                "recorder_dt":1,
                "n_phases":64,
                "hs_tau": 5.6,
                "hf_tau": 0.51,
                "conn_id":"1d",  
                "init_noise_seed":np.random.randint(0,100000),
                "record_handle_stell":{"stell_syn_inhib_g":{"state": True,"cells_to_record":"all"}},
                }
            mult_input_params[str(sim_num)] = input_params
            sim_num +=1

    return mult_input_params