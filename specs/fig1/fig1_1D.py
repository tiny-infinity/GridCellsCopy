import subprocess
import os
import time
import sys
import numpy as np



def generate_input_params():
    sim_dur = float(6000)
    input_params = {
        "sim_dur": sim_dur,
        "sim_id": "VALD-ACVT-1D-S-s-1a",
        "vel_type": "ACVT-1DAC",
        "intrnrn_noise":[sim_dur,0,2e-3],
        "stell_noise":[sim_dur,0,1e-3],
        "init_noise_seed":500,
        "noise_seed":1000,
        "n_cpus":40,
    }
    return input_params