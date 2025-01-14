import numpy as np




def generate_mult_input_params():
    sim_dur = float(80000)
    gh_range_arr = np.linspace(0,0.0015,16)
    dc_range_arr=np.array([ 1.60408103e-03,  1.46261512e-03,  1.28040971e-03,  1.10766164e-03,
        9.46567567e-04,  8.62474737e-04,  7.06476491e-04,  5.42020859e-04,
        3.73748911e-04,  2.07613748e-04,  3.95482557e-05, -1.32397520e-04,
       -2.99515027e-04, -4.69680325e-04, -6.41963358e-04, -8.21306155e-04])
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
