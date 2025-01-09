import numpy as np

def generate_input_params():
    sim_dur = float(12000)

    input_params = {        
        "sim_dur": sim_dur,
        "sim_id": "VALD-MPD-FSPK-S-s-1a",
        "traj_id": '201',
        "vel_type": "const",
        "g_h_bar":0.0015,
        "intrnrn_init_noise":[100,0,0.5],
        "stell_init_noise":[100,0,0.5],
        "stell_const_dc":[1.4e-3,-2.75e-3],
        "recorder_dt":0.025,
        "init_noise_seed":np.random.randint(100000),
        "record_handle_stell":{"stell_ih_mf":{"state": True,"cells_to_record":"all"},
                               "stell_ih_ms":{"state": True,"cells_to_record":"all"},
                               "stell_v":{"state": True,"cells_to_record":"all"},
                               "stell_syn_inhib_g":{"state": True,"cells_to_record":"all"}}

    }
    return input_params
