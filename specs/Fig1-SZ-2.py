import numpy as np

def generate_input_params():

    sim_dur = float(12000)
    mean_shift=24
    input_params = {
        "sim_dur": sim_dur,
        "sim_id": "VALD-ACVT-SIZE-S-s-2a",
        "ii_mean":mean_shift,
        "is_mean":mean_shift,
        "si_mean":mean_shift-10,

        "intrnrn_init_noise":[500,0,0.5],
        "stell_init_noise":[500,0,0.5],
        "intrnrn_noise":[sim_dur,0,2e-3],
        "stell_noise":[sim_dur,0,1e-3],
        "n_phases":int((mean_shift)*2),
        "init_noise_seed":500,
        "noise_seed":1000,
    }
    return input_params