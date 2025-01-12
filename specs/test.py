import subprocess
import os
import time
import sys
import numpy as np



def generate_input_params():
    sim_dur = float(12000)
    input_params = {
        "sim_dur": sim_dur,
        "sim_id": "test",
        "intrnrn_noise":[sim_dur,0,0],
        "stell_noise":[sim_dur,0,0],
        "init_noise_seed":500,
        # "stell_const_dc":[2e-3,-1],
        "stell_const_dc":[1.315e-3,-3e-2],
        "g_h_bar":0,
        # "noise_seed":1000,
        "n_cpus":40,
    }
    return input_params