

def generate_input_params():
    sim_dur = float(10000)

    input_params = {
        "sim_dur": sim_dur,
        "N_intrnrn": 1,
        "N_stell": 1,
        "N_per_sheet": 1,
        "sim_id": "VALD-ACVT-STDN-S-s-3a",
        "vel_type": "const",
        "is_peak":0,
        "intrnrn_init_noise":[100,0,0],
        "stell_init_noise":[100,0,0],
        "n_cpus":1,
        "intrnrn_dc_amp":1.5e-3,
        "stell_const_dc":[-2.66e-3,-1],
        "recorder_dt":0.025,
        "conn_id":"stdn",
        "instr_id":"stdn",
        "input_id":None,
        "omega_i_theta": 0.007,
        "init_noise_seed":50,
        "noise_seed":100,
        "record_handle_stell":{"stell_v":{"state": True,"cells_to_record":"all"}}
    }
    return input_params