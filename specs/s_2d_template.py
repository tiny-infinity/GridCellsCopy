import numpy as np


def generate_input_params():
    """
    Generates a dictionary of input parameters for a simulation.
    Returns:
        dict: A dictionary containing the provided input parameters:
    """
    sim_dur = 1000
    N_per_sheet = 2601
    N_per_axis=int(np.sqrt(N_per_sheet))
    N_sheets_stell=4
    N_stell = N_per_sheet*N_sheets_stell
    N_intrnrn = N_per_sheet
    
    input_params = {
        
        'sim_dur': float(sim_dur),
        'conn_id': "2d",
        "instr_id": "2d",
        'build_conn_matrix': False,
        'save_conn_matrix':True,
        'sim_id': 'BaseModel2D',
        "N_per_axis": N_per_axis,
        "N_per_sheet": N_per_sheet,
        "N_sheets_stell":N_sheets_stell,
        "N_stell": N_stell,
        "N_intrnrn": N_intrnrn,
        "is_peak": 2,
        "is_stdev": 1.05,
        "is_mean": 10,
        "si_peak": 1,
        "si_stdev": 1.05,
        "si_mean": [[7, 0],[-7,0],[0,7],[0,-7]],
        "stell_init_noise": [250,0,0.05],
        "intrnrn_init_noise": [250,0,0.05],
        "ii_peak": 2,
        "ii_stdev": 1.05,
        "ii_mean": 10,
        "stell_noise": [sim_dur, 0, 0.005],
        "intrnrn_dc_amp": 0.0005,
        "intrnrn_noise": [sim_dur, 0, 0.005],
        "stell_const_dc": [-4e-3,1e-4,1e-4,-4e-3]}
    return input_params
