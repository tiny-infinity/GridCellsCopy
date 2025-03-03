import numpy as np

def generate_mult_input_params():
    gh_range_arr= [0.0015,0]
    dc_range_arr = np.linspace(-2.7e-3,1e-2,500,endpoint=True)
    n_trials = 1
    sim_num = 0
    multiple_input_params = {}
    sim_dur = float(60000)
    for gh in gh_range_arr:
        for dc in dc_range_arr:
            for tr in range(n_trials):
                input_params = {
                    "sim_num":str(sim_num),
                    "sim_dur": sim_dur,
                    "sim_id": "VALD-HCN-SHRK-S-m-1b",
                    "intrnrn_init_noise":[100,0,0.5],
                    "stell_init_noise":[100,0,0.5],
                    "intrnrn_noise":[sim_dur,0,0],
                    "stell_noise":[sim_dur,0,0],
                    "stell_const_dc":[dc,-1],
                    "ii_peak": 1.4,
                    "g_h_bar":gh,
                    "init_noise_seed":np.random.randint(0,100000),
                    "noise_seed":np.random.randint(0,100000),
                    "n_nodes":10,
                    "data_root":"data/",
                }
                multiple_input_params[str(sim_num)] = input_params
                sim_num +=1
    return multiple_input_params