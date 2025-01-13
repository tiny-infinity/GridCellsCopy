
import numpy as np

def generate_mult_input_params():
    sim_dur = float(80000)
    dc_range_arr = np.linspace(-2.7e-3,1e-2,500,endpoint=False)
    gh_range_arr = np.linspace(0,0.0015,16)
    n_trials = 1
    sim_num = 0
    mult_input_params = {}
    for gh in gh_range_arr:
        for dc in dc_range_arr:
            for tr in range(n_trials):
                input_params = {
                    "sim_num":str(sim_num),
                    "sim_dur": sim_dur,
                    "sim_id": "VALD-PRED-INT-S-m-3a",
                    "vel_type": "PRED-IHD",
                    "g_h_bar":gh,
                    "data_root":"data/",
                    "n_nodes":12,
                    "intrnrn_init_noise":[100,0,0.5],
                    "stell_init_noise":[100,0,0.5],
                    "stell_const_dc":[1.4e-3,-2.75e-3],
                    "n_phases":64,
                    "allothetic_dur":3000,
                    "lambda0":2*np.pi,
                    "conn_id":"asym",   
                    "init_noise_seed":np.random.randint(0,100000),
                    "extra_params":{"stell_dc":dc,"dir_change_t":40000},
                    }
                mult_input_params[str(sim_num)] = input_params
                sim_num +=1

    return mult_input_params