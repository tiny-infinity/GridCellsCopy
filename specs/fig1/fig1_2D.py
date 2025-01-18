import numpy as np


def generate_input_params():
    sim_dur = 15000
    N_per_sheet = 5041
    N_per_axis=int(np.sqrt(N_per_sheet))
    N_sheets_stell=4
    N_stell = N_per_sheet*N_sheets_stell
    N_intrnrn = N_per_sheet
    si_shift = 6.25
    input_params = {
        'sim_dur': float(sim_dur),
        'conn_id': "2d",
        "instr_id": "2d",
        'build_conn_matrix': True,
        "save_conn_matrix":True,
        'sim_id': 'VALD-ACVT-2D-S-s-2a',
        "N_per_axis": N_per_axis,
        "N_per_sheet": N_per_sheet,
        "N_sheets_stell":N_sheets_stell,
        "N_stell": N_stell,
        "N_intrnrn": N_intrnrn,
        "is_stdev": 1.05,
        "is_mean": 10,
        "si_peak": 1,
        "is_peak": 2,
        "ii_peak": 1,
        "si_stdev": 1.05,
        "n_cpus":40,
        "si_mean": [[si_shift, 0],[-si_shift,0],[0,si_shift],[0,-si_shift]],
        "stell_init_noise": [1000,0,0.5],
        "intrnrn_init_noise": [1000,0,0.5],
        "noise_seed": 1000,
        "init_noise_seed": 5000,
        "ii_stdev": 1.05,
        "ii_mean": 10,
        "intrnrn_noise":[sim_dur,0,0],
        "stell_noise":[sim_dur,0,0],
        "intrnrn_dc_amp": 1e-3,
        "stell_const_dc": [0,-4e-3,-4e-3,-4e-3]}
    return input_params
