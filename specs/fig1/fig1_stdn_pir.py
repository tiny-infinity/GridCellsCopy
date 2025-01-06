import numpy as np

def generate_mult_input_params()-> dict:
    inhib_weight_arr = [0.15,0.175,0.25,0.65,2]
    n_trials=1
    sim_num = 0
    sim_dur = float(3500)
    mult_input_params = {}
    for i,inhib_weight in enumerate(inhib_weight_arr):
        for tr in range(n_trials):
            input_params = {
                "sim_dur": sim_dur,
                "sim_num": str(sim_num),
                "N_intrnrn": 1,
                "N_stell": 1,
                "N_per_sheet": 1,
                "sim_id": "VALD-ACVT-STDN-S-m-1a",
                "is_peak":0,
                "intrnrn_init_noise":[100,0,0],
                "stell_init_noise":[100,0,0],
                "n_cpus":1,
                "intrnrn_dc_amp":1.5e-3,
                "stell_const_dc":[-2.75e-3,-1],
                "recorder_dt":0.025,
                "conn_id":"stdn",
                "instr_id":"stdn",
                "input_id":"pir",
                "omega_i_theta": 0.007,
                "init_noise_seed":50,
                "noise_seed":100,
                "extra_params":{"start":1500,"weight":inhib_weight},
                "record_handle_stell":{"stell_v":{"state": True,"cells_to_record":"all"}}
            }
            
            mult_input_params[str(sim_num)] = input_params #Important.
            sim_num+=1
    return mult_input_params