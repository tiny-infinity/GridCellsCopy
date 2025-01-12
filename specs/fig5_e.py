import numpy as np


def generate_mult_input_params():
    sim_dur = float(300000)
    fast_tau_arr = np.linspace(0.51,2.5,10)
    slow_tau_arr = np.linspace(5.6,53.1,10)
    n_trials = 1
    sim_num = 0
    mult_input_params = {}
    for tau_s in slow_tau_arr:
        for tau_f in fast_tau_arr:
            for tr in range(n_trials):
                input_params = {
                    "sim_num":str(sim_num),
                    "sim_dur": sim_dur,
                    "sim_id": "VALD-PRED-INT-S-m-2b",
                    "vel_type": "PRED-IHD",
                    "g_h_bar":0.0015,
                    "n_nodes":10,
                    "data_root":"data/",
                    "netcon_delay":5,
                    "intrnrn_init_noise":[100,0,0.5],
                    "stell_init_noise":[100,0,0.5],
                    "stell_const_dc":[0.002,-2.75e-3],
                    "n_phases":64,
                    "allothetic_dur":3000,
                    "lambda0":2*np.pi,
                    "hs_tau":tau_s,
                    "hf_tau": tau_f,
                    "conn_id":"asym",   
                    "init_noise_seed":np.random.randint(0,100000),
                    "extra_params":{"stell_dc":0,"dir_change_t":150000}, # -0.0006799 0.0016496
                    }
                mult_input_params[str(sim_num)] = input_params
                sim_num +=1

    return mult_input_params
