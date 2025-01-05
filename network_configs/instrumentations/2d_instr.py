from neuron import h
from neuron.units import ms, mV
import numpy as np
import network_configs.instrumentations.instr_utils as instr_utils


    
def setup_instrumentation(network):
    params = network.params
    stell_l1_gids = list(range(0, params["N_per_sheet"]))
    stell_l2_gids = list(range(params["N_per_sheet"], 2*params["N_per_sheet"]))
    stell_l3_gids = list(range(2*params["N_per_sheet"], 3*params["N_per_sheet"]))
    stell_l4_gids = list(range(3*params["N_per_sheet"], 4*params["N_per_sheet"]))
    sim_dur = params["sim_dur"]
    for stell in network.stellate_cells:
        
        #DC input
        if stell._gid in stell_l1_gids:
            stell.ext_dc.dur = sim_dur
            stell.ext_dc.amp = network.params["stell_const_dc"][0]
        elif stell._gid in stell_l2_gids:
            stell.ext_dc.dur = sim_dur
            stell.ext_dc.amp = network.params["stell_const_dc"][1]
        elif stell._gid in stell_l3_gids:
            stell.ext_dc.dur = sim_dur
            stell.ext_dc.amp = network.params["stell_const_dc"][2]
        elif stell._gid in stell_l4_gids:
            stell.ext_dc.dur = sim_dur
            stell.ext_dc.amp = network.params["stell_const_dc"][3]
        
        instr_utils.set_stell_range_variables(stell,params)
        # Initial noise
        instr_utils.set_intial_noise(stell,params["stell_init_noise"])
        # Noise
        instr_utils.set_noise(stell,params["stell_noise"])
        #recorders
        instr_utils.setup_recorders(stell,params["record_handle_stell"],params["recorder_dt"])

    for interneuron in network.interneurons:
        interneuron.ext_dc.amp = params["intrnrn_dc_amp"]
        interneuron.ext_dc.dur = params["sim_dur"]
        
        instr_utils.set_intrnrn_range_variables(interneuron,params)        
        instr_utils.set_intial_noise(interneuron,params["intrnrn_init_noise"])
        instr_utils.set_noise(interneuron,params["intrnrn_noise"])
        instr_utils.setup_recorders(interneuron,params["record_handle_intrnrn"],params["recorder_dt"])
        




