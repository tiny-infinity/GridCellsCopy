from neuron import h
from neuron.units import ms, mV
import numpy as np
from trajectory import Trajectory
import sim_hf as s_hf
from network import Network


def setup_instrumentation(network):
    params = network.params

    #global variables
    h.hf_tau_input_stellate_mech = params['hf_tau']
    h.hs_tau_input_stellate_mech = params['hs_tau']
    h.phi_i_theta = params['phi_i_theta']  # 0
    h.omega_i_theta = params['omega_i_theta']  # 0.01
    
    stell_l1_gids = list(range(0,  network.N_per_sheet))
    stell_l2_gids = list(range( network.N_per_sheet,  network.N_per_sheet*2))
    network.traj = Trajectory(params)
    
    #vectors for DC input
    network.ext_t = h.Vector(np.arange(0, params["sim_dur"]+params["dt"], params["dt"]))
    network.ext_amp_right = h.Vector(network.traj.right_dc)
    network.ext_amp_left = h.Vector(network.traj.left_dc)
    network.ext_amp_intrnrn = h.Vector(network.traj.intrnrn_dc)

    for stell in network.stellate_cells:
        #range variables
        stell.soma(0.5).stellate_mech.ghbar = params['g_h_bar']
        stell.soma(0.5).stellate_mech.gnap_bar = params['gnap_bar']
        stell.soma(0.5).i_theta_stell.Amp = params['stell_theta_Amp']
        stell.soma(0.5).i_theta_stell.c = params['stell_theta_c']
        stell.soma(0.5).i_theta_stell.omega = params['stell_theta_omega']
        
        #DC input
        if stell._gid in stell_l1_gids:
            stell.ext_dc.dur = params['sim_dur']
            network.ext_amp_right.play(stell.ext_dc._ref_amp, network.ext_t, True)

        if stell._gid in stell_l2_gids:
            stell.ext_dc.dur = params['sim_dur']
            network.ext_amp_left.play(stell.ext_dc._ref_amp, network.ext_t, True)    

        #recorders
        s_hf.setup_recorders(stell,params["record_handle_stell"],params["recorder_dt"])


    for interneuron in network.interneurons:
        #range variables
        interneuron.soma(0.5).i_theta.Amp = params['Amp_i_theta']
        interneuron.ext_dc.dur = params['sim_dur'] - interneuron.ext_dc.delay

        #DC input
        network.ext_amp_intrnrn.play(interneuron.ext_dc._ref_amp, network.ext_t, True)
        
        #recorders
        s_hf.setup_recorders(interneuron,params["record_handle_intrnrn"],params["recorder_dt"])

