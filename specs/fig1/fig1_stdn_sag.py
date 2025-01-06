import numpy as np

def generate_mult_input_params()-> dict:
    hyper_dc_arr = [-3.5e-3,-5.5e-3,-7.5e-3,-9.5e-3,-1.15e-2]
    n_trials=1
    sim_num = 0
    sim_dur = float(1000)
    mult_input_params = {}
    for i,hyper_dc in enumerate(hyper_dc_arr):
        for tr in range(n_trials):
            input_params = {
                "sim_dur": sim_dur,
                "sim_num": sim_num,
                "N_intrnrn": 1,
                "N_stell": 1,
                "N_per_sheet": 1,
                "sim_id": "VALD-ACVT-STDN-S-m-2a",
                "is_peak":0,
                "intrnrn_init_noise":[100,0,0],
                "stell_init_noise":[100,0,0],
                "n_cpus":1,
                "intrnrn_dc_amp":1.5e-3,
                "stell_const_dc":[-2.65e-3,-1],
                "recorder_dt":0.025,
                "conn_id":"stdn",
                "instr_id":"stdn",
                "input_id":"sag",
                "omega_i_theta": 0.007,
                "init_noise_seed":50,
                "noise_seed":100,
                "extra_params":{"dur1":250,"dur2":1000,"amp1":-2.7e-3,"amp2":hyper_dc},
                "init_noise_seed":50,
                "record_handle_intrnrn":{ "intrnrn_v":{"state": True,"cells_to_record":"all"}},
                "record_handle_stell":{"stell_v":{"state": True,"cells_to_record":"all"}}
            }
            
            mult_input_params[str(sim_num)] = input_params #Important.
            sim_num+=1
    return mult_input_params