import subprocess
import os
import helper_functions as hf
import time
import numpy as np

"Finding a better DC value for zero velocity with stellate oscilations"


def def_sim_dict():
    m = 0.00571
    c = 0.0201
    stdev = 0.103
    mean = int((stdev-c)/m)
    stell_noise_std = np.linspace(0, 3*2e-3, 150)
    intrnrn_noise_std = np.linspace(0, 3*0.0005, 150)

    n_trials = 30
    sim_num = 0
    sim_dict = {}
    # for j in ss_width:
    for j in range(stell_noise_std.shape[0]):
        for tr in range(n_trials):
            input_params = {
                "sim_dur": float(30000),
                "sim_num": sim_num,
                "sim_id": "path_integ_955_noise",
                "intrnrn_dc_amp": 0.0005,
                "stell_noise": [30000, 0, stell_noise_std[j]],
                "stell_init_noise": [50, 0, 0.05],
                "intrnrn_noise": [30000, 0,  intrnrn_noise_std[j]],
                "intrnrn_init_noise": [50, 0, 0.05],
                "split_sim": [False, 1000],
                "vel_type": "input",
                "stell_theta_Amp": 0,
                "stell_theta_c": 0,
                "stell_theta_omega": 0.0276,
                "is_stdev": stdev,
                "ii_stdev": stdev,
                "allothetic_dur": 5000,
                "si_stdev": 1.21*stdev,
                "traj_id": 955,
                "ii_mean": mean,
                "is_mean": mean,
                "si_mean": 7}
            sim_dict[sim_num] = input_params
            sim_num += 1
    return sim_dict
