from neuron import h
from neuron.units import ms, mV
import numpy as np
from trajectory import Trajectory


    
def setup_instrumentation(network):
    params=network.params
    sim_dur = params['sim_dur']
    n_stell_pering = network.n_stell//2
    stell_l1_gids = list(range(0, n_stell_pering))
    stell_l2_gids = list(range(n_stell_pering, n_stell_pering*2))
    traj = Trajectory(params)
    dc_vel_tuning=0 #np.max(np.concatenate((traj.allothetic_right_dc,traj.allothetic_left_dc)))
    vel_type = params['vel_type']
    network.ext_t = h.Vector(np.arange(0, sim_dur+h.dt, h.dt))
    network.ext_amp_right = h.Vector(traj.right_dc)
    network.ext_amp_left = h.Vector(traj.left_dc)
    network.ext_amp_intrnrn = h.Vector(traj.intrnrn_dc)
    recorder_dt = params["recorder_dt"]
    # print(list(network.ext_amp_intrnrn))
    if vel_type == 'input':
        network.ext_amp_right_allo = h.Vector(traj.allothetic_right_dc)
        network.ext_amp_left_allo = h.Vector(traj.allothetic_left_dc)
        network.ext_amp_intnrn_allo = h.Vector(traj.allothetic_intrnrn_dc)


    
        

    for stell in network.stellate_cells:
        #electrodes
        stell.soma(0.5).stellate_mech.ghbar = params['g_h_bar']
        stell.soma(0.5).stellate_mech.gnap_bar = params['gnap_bar']

        stell.soma(0.5).i_theta_stell.Amp = params['stell_theta_Amp']
        stell.soma(0.5).i_theta_stell.c = params['stell_theta_c']
        stell.soma(
            0.5).i_theta_stell.omega = params['stell_theta_omega']
        
        if stell._gid in stell_l1_gids:
            stell.ext_dc.dur = sim_dur
            if stell._gid in traj.active_cells and vel_type == 'input':
                network.ext_amp_right_allo.play(stell.ext_dc._ref_amp, network.ext_t, True)
            else:
                
                network.ext_amp_right.play(stell.ext_dc._ref_amp, network.ext_t, True)
                

        if stell._gid in stell_l2_gids:
            stell.ext_dc.dur = sim_dur
            # allothetic input
            if stell._gid in traj.active_cells and vel_type == 'input':
                network.ext_amp_left_allo.play(stell.ext_dc._ref_amp, network.ext_t, True)
            else:
                network.ext_amp_left.play(stell.ext_dc._ref_amp, network.ext_t, True)    

        #recorders
        if params['record_handle_stell']['stell_v']['state']==True:
            if params['record_handle_stell']['stell_v']['cells_to_record']=='all':
                stell.recorder['stell_v'] = h.Vector().record(stell.soma(0.5)._ref_v,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_v']['cells_to_record']):
                stell.recorder['stell_v'] = h.Vector().record(stell.soma(0.5)._ref_v,recorder_dt)
        
        if params['record_handle_stell']['stell_spks']['state']==True:
            if params['record_handle_stell']['stell_spks']['cells_to_record']=='all':
                stell.recorder['stell_spks'] = h.Vector()  # spikes
                stell._spike_detector.record(stell.recorder['stell_spks'])
            elif stell._gid in list(params['record_handle_stell']['stell_spks']['cells_to_record']):
                stell.recorder['stell_spks'] = h.Vector()  # spikes
                stell._spike_detector.record(stell.recorder['stell_spks'])

        if params['record_handle_stell']['stell_ih']['state']==True:
            if params['record_handle_stell']['stell_ih']['cells_to_record']=='all':
                stell.recorder['stell_ih'] = h.Vector().record(stell.soma(0.5).stellate_mech._ref_ih,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_ih']['cells_to_record']):
                stell.recorder['stell_ih'] = h.Vector().record(stell.soma(0.5).stellate_mech._ref_ih,recorder_dt)

        if params['record_handle_stell']['stell_ina']['state']==True:
            if params['record_handle_stell']['stell_ina']['cells_to_record']=='all':
                stell.recorder['stell_ina'] = h.Vector().record(stell.soma(0.5)._ref_ina,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_ina']['cells_to_record']):
                stell.recorder['stell_ina'] = h.Vector().record(stell.soma(0.5)._ref_ina,recorder_dt)

        if params['record_handle_stell']['stell_ih_mf']['state']==True:
            if params['record_handle_stell']['stell_ih_mf']['cells_to_record']=='all':
                stell.recorder['stell_ih_mf'] = h.Vector().record(stell.soma(0.5).stellate_mech._ref_mhf,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_ih_mf']['cells_to_record']):
                stell.recorder['stell_ih_mf'] = h.Vector().record(stell.soma(0.5).stellate_mech._ref_mhf,recorder_dt)

        if params['record_handle_stell']['stell_ih_ms']['state']==True:
            if params['record_handle_stell']['stell_ih_ms']['cells_to_record']=='all':
                stell.recorder['stell_ih_ms'] = h.Vector().record(stell.soma(0.5).stellate_mech._ref_mhs,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_ih_ms']['cells_to_record']):
                stell.recorder['stell_ih_ms'] = h.Vector().record(stell.soma(0.5).stellate_mech._ref_mhs,recorder_dt)

        if params['record_handle_stell']['stell_ik']['state']==True:
            if params['record_handle_stell']['stell_ik']['cells_to_record']=='all':
                stell.recorder['stell_ik'] = h.Vector().record(stell.soma(0.5)._ref_ik,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_ik']['cells_to_record']):
                stell.recorder['stell_ik'] = h.Vector().record(stell.soma(0.5)._ref_ik,recorder_dt)
        if params['record_handle_stell']['stell_gna']['state']==True:
            if params['record_handle_stell']['stell_gna']['cells_to_record']=='all':
                stell.recorder['stell_gna'] = h.Vector().record(stell.soma(0.5).stellate_mech._ref_gna,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_gna']['cells_to_record']):
                stell.recorder['stell_gna'] = h.Vector().record(stell.soma(0.5).stellate_mech._ref_gna,recorder_dt)

        if params['record_handle_stell']['stell_syn_inhib_i']['state']==True:
            if params['record_handle_stell']['stell_syn_inhib_i']['cells_to_record']=='all':
                stell.recorder['stell_syn_inhib_i'] = h.Vector().record(stell.inhb_syn._ref_i,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_syn_inhib_i']['cells_to_record']):
                stell.recorder['stell_syn_inhib_i'] = h.Vector().record(stell.inhb_syn._ref_i,recorder_dt)

        if params['record_handle_stell']['stell_syn_inhib_g']['state']==True:
            if params['record_handle_stell']['stell_syn_inhib_g']['cells_to_record']=='all':
                stell.recorder['stell_syn_inhib_g'] = h.Vector().record(stell.inhb_syn._ref_g,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_syn_inhib_g']['cells_to_record']):
                stell.recorder['stell_syn_inhib_g'] = h.Vector().record(stell.inhb_syn._ref_g,recorder_dt)

        if params['record_handle_stell']['stell_ext_dc_i']['state']==True:
            if params['record_handle_stell']['stell_ext_dc_i']['cells_to_record']=='all':
                stell.recorder['stell_ext_dc_i'] = h.Vector().record(stell.ext_dc._ref_i,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_ext_dc_i']['cells_to_record']):
                stell.recorder['stell_ext_dc_i'] = h.Vector().record(stell.ext_dc._ref_i,recorder_dt)

        if params['record_handle_stell']['stell_ext_dc_amp']['state']==True:
            if params['record_handle_stell']['stell_ext_dc_amp']['cells_to_record']=='all':
                stell.recorder['stell_ext_dc_amp'] = h.Vector().record(stell.ext_dc._ref_amp,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_ext_dc_amp']['cells_to_record']):
                stell.recorder['stell_ext_dc_amp'] = h.Vector().record(stell.ext_dc._ref_amp,recorder_dt)
        
        if params['record_handle_stell']['stell_theta_i']['state']==True:
            if params['record_handle_stell']['stell_theta_i']['cells_to_record']=='all':
                stell.recorder['stell_theta_i'] = h.Vector().record(stell.soma(0.5).i_theta_stell._ref_itheta,recorder_dt)
            elif stell._gid in list(params['record_handle_stell']['stell_theta_i']['cells_to_record']):
                stell.recorder['stell_theta_i'] = h.Vector().record(stell.soma(0.5).i_theta_stell._ref_itheta,recorder_dt)

    for interneuron in network.interneurons:
        #electrodes
        # interneuron.soma(0.5).i_theta.dc=max(params['stell_const_dc'])
        # interneuron.ext_dc.amp = params['intrnrn_dc_amp']
        if params['vel_type'] == 'input':
            # interneuron.ext_dc.delay =params["intrnrn_init_noise"][0]
            # interneuron.soma(0.5).i_theta.dc = dc_vel_tuning
            pass

        else:
            # interneuron.ext_dc.delay = 0
            # interneuron.soma(0.5).i_theta.dc = np.max(params['stell_const_dc'])
            pass

        interneuron.ext_dc.dur = params['sim_dur'] - \
            interneuron.ext_dc.delay
        if interneuron._gid in traj.active_cells and vel_type == 'input':
            network.ext_amp_intnrn_allo.play(interneuron.ext_dc._ref_amp, network.ext_t, True)
        
        elif interneuron._gid==params["n_intrnrn"]+params["n_stell"]-1:
            network.intrnrn_const_dc = h.Vector(np.full_like(network.ext_t,params["global_inhib_dc"]))
            network.intrnrn_const_dc.play(interneuron.ext_dc._ref_amp, network.ext_t, True)
        else:
            network.ext_amp_intrnrn.play(interneuron.ext_dc._ref_amp, network.ext_t, True)
            
            
        
        interneuron.soma(0.5).i_theta.Amp = params['Amp_i_theta']
        
        #recorders
        if params['record_handle_intrnrn']['intrnrn_spks']['state']==True:
            if params['record_handle_intrnrn']['intrnrn_spks']['cells_to_record']=='all':
                interneuron.recorder['intrnrn_spks'] = h.Vector()  # spikes
                interneuron._spike_detector.record(interneuron.recorder['intrnrn_spks'])
            elif interneuron._gid in list(params['record_handle_intrnrn']['intrnrn_spks']['cells_to_record']):
                interneuron.recorder['intrnrn_spks']= h.Vector()  # spikes
                interneuron._spike_detector.record(interneuron.recorder['intrnrn_spks'])

        if params['record_handle_intrnrn']['intrnrn_v']['state']==True:
            if params['record_handle_intrnrn']['intrnrn_v']['cells_to_record']=='all':
                interneuron.recorder['intrnrn_v'] = h.Vector().record(interneuron.soma(0.5)._ref_v,recorder_dt)
            elif interneuron._gid in list(params['record_handle_intrnrn']['intrnrn_v']['cells_to_record']):
                interneuron.recorder['intrnrn_v'] = h.Vector().record(interneuron.soma(0.5)._ref_v,recorder_dt)

        if params['record_handle_intrnrn']['intrnrn_syn_exc_g']['state']==True:
            if params['record_handle_intrnrn']['intrnrn_syn_exc_g']['cells_to_record']=='all':
                interneuron.recorder['intrnrn_syn_exc_g'] = h.Vector().record(interneuron.exc_syn._ref_g,recorder_dt)
            elif interneuron._gid in list(params['record_handle_intrnrn']['intrnrn_syn_exc_g']['cells_to_record']):
                interneuron.recorder['intrnrn_syn_exc_g'] = h.Vector().record(interneuron.exc_syn._ref_g,recorder_dt)

        if params['record_handle_intrnrn']['intrnrn_syn_inhib_g']['state']==True:
            if params['record_handle_intrnrn']['intrnrn_syn_inhib_g']['cells_to_record']=='all':
                interneuron.recorder['intrnrn_syn_inhib_g'] = h.Vector().record(interneuron.inhb_syn._ref_g,recorder_dt)
            elif interneuron._gid in list(params['record_handle_intrnrn']['intrnrn_syn_inhib_g']['cells_to_record']):
                interneuron.recorder['intrnrn_syn_inhib_g'] = h.Vector().record(interneuron.inhb_syn._ref_g,recorder_dt)
        
        if params['record_handle_intrnrn']['intrnrn_ext_dc_i']['state']==True:
            if params['record_handle_intrnrn']['intrnrn_ext_dc_i']['cells_to_record']=='all':
                interneuron.recorder['intrnrn_ext_dc_i'] = h.Vector().record(interneuron.ext_dc._ref_i,recorder_dt)
            elif interneuron._gid in list(params['record_handle_intrnrn']['intrnrn_ext_dc_i']['cells_to_record']):
                interneuron.recorder['intrnrn_ext_dc_i'] = h.Vector().record(interneuron.ext_dc._ref_i,recorder_dt)

        if params['record_handle_intrnrn']['intrnrn_theta_i']['state']==True:
            if params['record_handle_intrnrn']['intrnrn_theta_i']['cells_to_record']=='all':
                interneuron.recorder['intrnrn_theta_i'] = h.Vector().record(interneuron.soma(0.5).i_theta._ref_itheta,recorder_dt)
            elif interneuron._gid in list(params['record_handle_intrnrn']['intrnrn_theta_i']['cells_to_record']):
                interneuron.recorder['intrnrn_theta_i'] = h.Vector().record(interneuron.soma(0.5).i_theta._ref_itheta,recorder_dt)
        
        if params['record_handle_intrnrn']['intrnrn_ext_dc_amp']['state']==True:
            if params['record_handle_intrnrn']['intrnrn_ext_dc_amp']['cells_to_record']=='all':
                interneuron.recorder['intrnrn_ext_dc_amp'] = h.Vector().record(interneuron.ext_dc._ref_amp,recorder_dt)
            elif interneuron._gid in list(params['record_handle_intrnrn']['intrnrn_ext_dc_amp']['cells_to_record']):
                interneuron.recorder['intrnrn_ext_dc_amp'] = h.Vector().record(interneuron.ext_dc._ref_amp,recorder_dt)



