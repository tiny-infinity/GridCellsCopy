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
            "sim_id": "VALD-VI-TRAJ-S-m-1a", #"VI-TI-SIMS-S-1a"
            "traj_id": '502',
            "vel_type": "input",
            "init_allothetic_input": True,
            "Amp_i_theta":8e-5,
            "allothetic_stell_dc":-0.0027,
            "intrnrn_init_noise":[100,0,0.5],
            "stell_init_noise":[100,0,0.5],
            "intrnrn_noise":[sim_dur,0,2e-3],
            "stell_noise":[sim_dur,0,1e-3],
            "stell_const_dc":[-2.453e-3,-2.75e-3],
            "n_phases":64,
            "vel_integ_zero":-0.002633,
            "vel_integ_or":-0.002906,
            "lambda0":2*np.pi,
            "allothetic_nrn_n":10,
            "init_noise_seed":np.random.randint(0,100000),
            "noise_seed":np.random.randint(0,100000),
            "n_cpus":40,
            "record_handle_stell":{"stell_syn_inhib_g":{"state":True,"cells_to_record":[46,64]},
                                     "stell_ext_dc_amp":{"state":True,"cells_to_record":[46,64]}}
        }
        multiple_input_params[str(sim_num)] = input_params
        sim_num +=1
    return multiple_input_params
