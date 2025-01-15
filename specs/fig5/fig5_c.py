import numpy as np




def generate_mult_input_params():
    sim_dur = float(80000)
    gh_range_arr = np.linspace(0,0.0015,16)
    dc_range_arr=np.array([ 1.28031834e-03,  1.04879480e-03,  8.16224066e-04,  5.64521873e-04,
        3.11440625e-04,  5.91552781e-05, -2.10812800e-04, -4.46674041e-04,
       -6.97370230e-04, -9.49550210e-04, -1.17218007e-03, -1.41407774e-03,
       -1.64949788e-03, -1.87639332e-03, -2.09076233e-03, -2.30484610e-03])
    n_trials = 1
    sim_num = 0
    mult_input_params = {}
    for i,gh in enumerate(gh_range_arr):
        for tr in range(n_trials):
            input_params = {
                "sim_num":str(sim_num),
                "sim_dur": sim_dur,
                "n_nodes":10,
                "data_root":"data/",
                "sim_id": "VALD-PRED-INT-S-m-4a",
                "vel_type": "PRED-IHD",
                "g_h_bar":gh,
                "intrnrn_init_noise":[100,0,0.5],
                "stell_init_noise":[100,0,0.5],
                "intrnrn_dc_amp":1.5e-3,
                "stell_const_dc":[0.002,-2.75e-3],
                "n_phases":64,
                "allothetic_dur":3000,
                "lambda0":2*np.pi,
                "hs_tau":5.6,
                "hf_tau": 0.51,
                "conn_id":"asym",   
                "init_noise_seed":np.random.randint(0,100000),
                "extra_params":{"stell_dc":dc_range_arr[i],"dir_change_t":40000}, # -0.0006799 0.0016496
                }
            mult_input_params[str(sim_num)] = input_params
            sim_num +=1

    return mult_input_params
