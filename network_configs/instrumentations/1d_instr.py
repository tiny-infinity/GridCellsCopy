"""
Instrumentation for 1D network
1. Set global variables and range variables.
2. Set initial and simulation noise
3. Setup recorders.
4. Setup external inputs by intializing Trajectory1D and playing it into the 
    network IClamps.
"""
from neuron import h
from neuron.units import ms, mV
import numpy as np
from network_configs.instrumentations.trajectory1D import Trajectory1D
import network_configs.instrumentations.instr_utils as  instr_utils


def setup_instrumentation(network):
    params = network.params
    
    instr_utils.set_global_variables(params)
    stell_l1_gids = list(range(0,  params["N_per_sheet"]))
    stell_l2_gids = list(range(params["N_per_sheet"],params["N_per_sheet"]*2))
    network.traj = Trajectory1D(params)
    """
        # --- Cue-based modulation of theta oscillation amplitude ---
    pos_array = np.array(network.traj.pos_rinb)           # position trace
    dc_input  = np.array(network.traj.intrnrn_dc)         # full input trace from trajectory
    
    # Parameters from params
    cue_pos  = params.get("cue_pos", 100.0)    # cm
    sigma    = params.get("cue_sigma", 15.0)   # cm
    L        = params.get("track_length", np.max(pos_array))
    
    # Helper: circular distance
    def periodic_distance(x, y, L):
        return np.minimum(np.abs(x - y), L - np.abs(x - y))
    
    # Compute Gaussian bump modulation
    dist = periodic_distance(pos_array, cue_pos, L)
    amp_mod = np.exp(-(dist**2) / (2 * sigma**2))   # [0,1]
    
    # Separate baseline + oscillation
    baseline = np.mean(dc_input)                    # constant offset
    oscill   = dc_input - baseline                  # zero-mean oscillatory component
    
    # Apply modulation only to oscillatory component
    modulated_intrnrn_dc = baseline + oscill * (1 - amp_mod)

    
    # Replace vector (this will be played into interneurons below)
    network.ext_amp_intrnrn = h.Vector(modulated_intrnrn_dc)
    """
    #create common vectors
    network.ext_t = h.Vector(np.arange(0, params["sim_dur"]+params["dt"], params["dt"]))
    network.ext_amp_right = h.Vector(network.traj.right_dc)
    network.ext_amp_left = h.Vector(network.traj.left_dc)
    network.ext_amp_intrnrn = h.Vector(network.traj.intrnrn_dc)

    for stell in network.stellate_cells:
        #range variables
        instr_utils.set_stell_range_variables(stell,params)

        if stell._gid in stell_l1_gids:
            stell.ext_dc.dur = params["sim_dur"]
            network.ext_amp_right.play(stell.ext_dc._ref_amp, network.ext_t, True)
                
        if stell._gid in stell_l2_gids:
            stell.ext_dc.dur = params["sim_dur"]
            #init allothetic input
            network.ext_amp_left.play(stell.ext_dc._ref_amp, network.ext_t, True)

        
        # Initial noise
        instr_utils.set_intial_noise(stell,params["stell_init_noise"],noise_seed=params["init_noise_seed"])
        # Noise
        instr_utils.set_noise(stell,params["stell_noise"],noise_seed=params["noise_seed"])
        #recorders
        instr_utils.setup_recorders(stell,params["record_handle_stell"],params["recorder_dt"])



    for interneuron in network.interneurons:
        #range variables
        interneuron.ext_dc.dur = params['sim_dur'] - interneuron.ext_dc.delay

        if interneuron._gid in network.traj.active_cells and params["vel_type"] == 'input':
            network.traj.ext_amp_intnrn_allo.play(interneuron.ext_dc._ref_amp, network.ext_t, True)
        else:
            network.ext_amp_intrnrn.play(interneuron.ext_dc._ref_amp, network.ext_t, True)
        instr_utils.set_intrnrn_range_variables(interneuron,params)
        instr_utils.set_intial_noise(interneuron,params["intrnrn_init_noise"],noise_seed=params["init_noise_seed"])
        instr_utils.set_noise(interneuron,params["intrnrn_noise"],noise_seed=params["noise_seed"])
        instr_utils.setup_recorders(interneuron,params["record_handle_intrnrn"],params["recorder_dt"])

