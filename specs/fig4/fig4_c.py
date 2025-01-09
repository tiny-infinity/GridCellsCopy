import numpy as np

def generate_mult_input_params():
    dc_range_arr = np.linspace(-2.7e-3,1e-2,1000,endpoint=True)
    coeffs=np.array([-2.56290623e+03,  6.04362519e+01, -4.62960821e-01,  3.14888925e-03])
    n_trials = 1
    sim_num = 0
    multiple_input_params = {}
    sim_dur = float(40000)
    for i,dc in enumerate(dc_range_arr):
        for tr in range(n_trials):
            input_params = {
                "sim_num":str(sim_num),
                "sim_dur": sim_dur,
                "sim_id": "VALD-HCN-SHRK-S-m-2a",
                "intrnrn_init_noise":[100,0,0.5],
                "stell_init_noise":[100,0,0.5],
                "intrnrn_noise":[sim_dur,0,0],
                "stell_noise":[sim_dur,0,0],
                "stell_const_dc":[dc,-1],
                "g_h_bar":0,
                "intrnrn_dc_amp":np.polyval(coeffs,dc),
                "init_noise_seed":np.random.randint(0,100000),
                "noise_seed":np.random.randint(0,100000),
                "n_nodes":10,
                "data_root":"data/",
            }
            multiple_input_params[str(sim_num)] = input_params
            sim_num +=1
    return multiple_input_params