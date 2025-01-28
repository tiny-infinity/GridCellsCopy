"""
Instrumentation for single cell simulations for Fig1A
"""
from neuron import h
from neuron.units import ms, mV
import numpy as np
import network_configs.instrumentations.instr_utils as  instr_utils

def create_piecewise(x,y,sim_dur):
    '''input x: time at step changes (including endpoints),y: step values'''
    t_ = np.arange(0,sim_dur+h.dt,h.dt)
    l = []
    for i in range(len(x)-1):
        l.append(np.logical_and(t_>=x[i],t_<=x[i+1]))
        
    return np.piecewise(t_,l,y)
    
def setup_instrumentation(network):
    params=network.params
    sim_dur = params['sim_dur']
    # instr_utils.set_global_variables(params)
    network.ext_t = h.Vector(np.arange(0, sim_dur+h.dt, h.dt))
    recorder_dt = params["recorder_dt"]
    network.ext_amp_stell = h.Vector(np.full_like( network.ext_t, network.params["stell_const_dc"][0]))
    network.ext_amp_intrnrn = h.Vector(np.full_like( network.ext_t, network.params["intrnrn_dc_amp"]))
    
    if network.params["input_id"]:
        add_xtra_inputs(network)
    
    for stell in network.stellate_cells:
        instr_utils.set_stell_range_variables(stell,params)
        stell.ext_dc.dur = sim_dur
        network.ext_amp_stell.play(stell.ext_dc._ref_amp, network.ext_t, True)

        #recorders
        instr_utils.setup_recorders(stell,params["record_handle_stell"],recorder_dt)

    for interneuron in network.interneurons:
        interneuron.ext_dc.dur = sim_dur
        network.ext_amp_intrnrn.play(interneuron.ext_dc._ref_amp, network.ext_t, True)
        instr_utils.set_intrnrn_range_variables(interneuron,params)
        
        #recorders
        instr_utils.setup_recorders(interneuron,params["record_handle_intrnrn"],recorder_dt)




def input_profiles(params):
    if params["input_id"]=="sag":
        dur1 = params["extra_params"]["dur1"]
        dur2 = params["extra_params"]["dur2"]
        amp1 = params["extra_params"]["amp1"]
        amp2 = params["extra_params"]["amp2"]
        x = [0,dur1,dur2]
        y= [amp1,amp2]
        return create_piecewise(x,y,params["sim_dur"])

def add_net_stim(network):
    network.ns = h.NetStim()
    network.ns.number = 1
    network.ns.start = network.params["extra_params"]["start"]
    network.nc = h.NetCon(network.ns, network.stellate_cells[0].inhb_syn)
    network.nc.weight[0]=network.params["extra_params"]["weight"]*1e-4

def add_xtra_inputs(network):
    if network.params["input_id"]=="sag":
        network.input_amp_stell = h.Vector(input_profiles(network.params))
        for stell in network.stellate_cells:
            stell.instr["IClamps"].append(h.IClamp(stell.soma(0.5)))
            stell.instr["IClamps"][-1].dur = network.params["sim_dur"]
            network.input_amp_stell.play(stell.instr["IClamps"][-1]._ref_amp, network.ext_t, True)
    if network.params["input_id"]=="pir":
        add_net_stim(network)
    if network.params["input_id"]=="resonance":
        sim_dur=network.params["sim_dur"]
        from scipy.signal import chirp
        delay = 2000
        t = np.linspace(0, sim_dur/1000, int((sim_dur-delay)/0.025)+1)
        network.vec_t=h.Vector(np.arange(delay, sim_dur+h.dt, h.dt))
        network.chirp_input = 1e-6*chirp(t, f0=network.params["extra_params"]["f0"], f1=network.params["extra_params"]["f1"], t1=sim_dur/1000, method='linear')
        network.chirp_input=h.Vector(network.chirp_input)
        for stell in network.stellate_cells:
            stell.instr["IClamps"].append(h.IClamp(stell.soma(0.5)))
            stell.instr["IClamps"][-1].dur = sim_dur
            stell.instr["IClamps"][-1].delay = delay

            network.chirp_input.play(stell.instr["IClamps"][-1]._ref_amp,network.vec_t, True)




