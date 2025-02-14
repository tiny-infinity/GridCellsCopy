"""
Template specs file to run 1D simulation.
"""
import numpy as np

def generate_input_params()-> dict:
    """Generates a dictionary of input parameters for a simulation.
    
    """
    sim_dur=8000 #in ms
    input_params = {   
        "sim_dur": sim_dur,
        "N_stell":2,
        "N_intrnrn":2,   
        "sim_id": "MotifSwitch",
        "conn_id":"motif",
        "instr_id":"motif",
        "input_id":"pulse",
        "extra_params":{"pulse_width":20,"amp":1e-2,"start":1000,"ipi":2000,
                        "first_cell_input":0},
        "init_noise_seed":500,
        "is_peak":1,
        "si_peak":1,
        "ii_peak":1,
        "stell_const_dc":-2.45e-3,
        "intrnrn_dc_amp":1e-3,
        "show_progress_bar":True,
        "n_cpus":4, #all avail
        "split_sim":[False,1000],
        "recorder_dt":0.025,
        "record_handle_stell": {
            "stell_v": {
                "state": True,
                "cells_to_record": "all"
            }},
        "record_handle_intrnrn": {
            "intrnrn_v": {
                "state": True,
                "cells_to_record": "all"
            },
            "intrnrn_switch_pulse": {
            "state": True,
            "cells_to_record": "all",
        }},
    }
    return input_params
