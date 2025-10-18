import numpy as np

def generate_mult_input_params():
    n_trials = 10
    sim_num = 0
    multiple_input_params = {}
    sim_dur = float(60000)
    for tr in range(n_trials):
        input_params = {
            "N_intrnrn": 192,
            "N_stell": 384,
            "N_per_sheet": 192,
            "is_stdev":0.103,
            "is_mean":14,
            "si_peak": 0.5,
            "si_stdev": 0.12462999999999999,
            "ii_stdev": 0.103,
            "ii_mean": 14,
            "sim_num":str(sim_num),
            "sim_dur": sim_dur,
            "sim_id": "NEW_THETA_1",
            "traj_id": 'tmod2',
            "stell_theta_omega": 0.0276,
            "vel_type": "input",
            "init_allothetic_input": True,
            "allothetic_stell_dc":-0.0027,
            "intrnrn_init_noise":[100,0,0.05],
            "stell_const_dc": [-1,-0.002],
            "phi_i_theta": 1.57,
            "hs_tau": 0,
            "hf_tau": 0,
            "stell_init_noise":[100,0,0.5],
            "intrnrn_noise":[sim_dur,0,2e-3],
            "stell_noise":[sim_dur,0,1e-3],
            "stell_const_dc":[-2.453e-3,-2.75e-3],
            "n_phases":64,
            "vel_integ_or":-0.002906,
            "lambda0":2*np.pi,
            "allothetic_nrn_n":10,
            "Amp_i_theta":0,
            "intrnrn_dc_amp":1e-3,
            "init_noise_seed":np.random.randint(0,100000),
            "noise_seed":np.random.randint(0,100000),
            "n_cpus":4,
            "tuning":0,
            "netcon_delay":1,
            "record_handle_stell":{"stell_syn_inhib_g":{"state":False,"cells_to_record":[46,64]},
                                     "stell_ext_dc_amp":{"state":False,"cells_to_record":[46,64]}},
            "record_handle_intrnrn":{"intrnrn_v":{"state": True,"cells_to_record":"all"}}
        }
        multiple_input_params[str(sim_num)] = input_params
        sim_num +=1
    return multiple_input_params
