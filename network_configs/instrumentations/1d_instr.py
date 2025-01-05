from neuron import h
from neuron.units import ms, mV
import numpy as np
from network_configs.instrumentations.trajectory1D import Trajectory1D
import network_configs.instrumentations.instr_utils as  instr_utils
from network import Network


def setup_instrumentation(network):
    params = network.params
    
    instr_utils.set_global_variables(params)
    stell_l1_gids = list(range(0,  params["N_per_sheet"]))
    stell_l2_gids = list(range(params["N_per_sheet"],params["N_per_sheet"]*2))
    network.traj = Trajectory1D(params)
    
    #create common vectors
    network.ext_t = h.Vector(np.arange(0, params["sim_dur"]+params["dt"], params["dt"]))
    network.ext_amp_right = h.Vector(network.traj.right_dc)
    network.ext_amp_left = h.Vector(network.traj.left_dc)
    network.ext_amp_intrnrn = h.Vector(network.traj.intrnrn_dc)

    for stell in network.stellate_cells:
        #range variables
        instr_utils.set_stell_range_variables(stell,params)
        
        #DC input
        if stell._gid in stell_l1_gids:
            stell.ext_dc.dur = params['sim_dur']
            network.ext_amp_right.play(stell.ext_dc._ref_amp, network.ext_t, True)

        if stell._gid in stell_l2_gids:
            stell.ext_dc.dur = params['sim_dur']
            network.ext_amp_left.play(stell.ext_dc._ref_amp, network.ext_t, True)    
        
        # Initial noise
        instr_utils.set_intial_noise(stell,params["stell_init_noise"])
        # Noise
        instr_utils.set_noise(stell,params["stell_noise"])
        #recorders
        instr_utils.setup_recorders(stell,params["record_handle_stell"],params["recorder_dt"])



    for interneuron in network.interneurons:
        #range variables
        interneuron.ext_dc.dur = params['sim_dur'] - interneuron.ext_dc.delay
        network.ext_amp_intrnrn.play(interneuron.ext_dc._ref_amp, network.ext_t, True)
        
        instr_utils.set_intrnrn_range_variables(interneuron,params)
        instr_utils.set_intial_noise(interneuron,params["intrnrn_init_noise"])
        instr_utils.set_noise(interneuron,params["intrnrn_noise"])
        instr_utils.setup_recorders(interneuron,params["record_handle_intrnrn"],params["recorder_dt"])

