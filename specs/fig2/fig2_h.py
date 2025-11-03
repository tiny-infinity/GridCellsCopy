import numpy as np

def generate_mult_input_params():
    n_trials = 1
    sim_num = 0
    multiple_input_params = {}
    sim_dur = float(30000)
    for tr in range(n_trials):
        input_params = {
            "sim_num":str(sim_num),
            "sim_dur": sim_dur,
            "sim_id": "fig2_h",
            "traj_id": '1',
            "vel_type": "input",
            "init_allothetic_input": True,
            "allothetic_stell_dc":-0.0027,
            "intrnrn_init_noise":[100,0,0.5],
            "stell_init_noise":[100,0,0.5],
            "intrnrn_noise":[sim_dur,0,2e-3],
            "stell_noise":[sim_dur,0,1e-3],
            "stell_const_dc":[-2.453e-3,-2.75e-3],
            "dc_amp_baseline":0,
            "n_phases":64,
            "vel_integ_or":-0.002906,
            "lambda0":2*np.pi,
            "allothetic_nrn_n":10,
            "Amp_i_theta":0,
            "intrnrn_dc_amp":0,
            "init_noise_seed":np.random.randint(0,100000),
            "noise_seed":np.random.randint(0,100000),
            "n_cpus":4,
            "netcon_delay":1,
            "tuning":0,
            #"record_handle_stell":{"stell_syn_inhib_g":{"state":False,"cells_to_record":[46,64]},
                                    # "stell_ext_dc_amp":{"state":False,"cells_to_record":[46,64]}},
            "record_handle_intrnrn":{"intrnrn_theta_dc":{"state":True,
        "cells_to_record":"all",
        "loc":"ext_theta_dc:_ref_amp"
    },
                                     "intrnrn_v": {"state": True,"cells_to_record": "all","loc":"soma(0.5):_ref_v"},
                                     "intrnrn_ext_dc_i": {"state": True,"cells_to_record": "all","loc":"ext_dc:_ref_i"},
                                     "intrnrn_ext_dc_amp": {"state": True,"cells_to_record":"all","loc":"ext_dc:_ref_amp"},
                                     "intrnrn_theta_i": {"state": True,"cells_to_record":"all","loc":"soma(0.5):i_theta:_ref_itheta"}}
        }
        multiple_input_params[str(sim_num)] = input_params
        sim_num +=1
    return multiple_input_params
