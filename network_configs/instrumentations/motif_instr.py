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
    network.ext_amp_stell = h.Vector(np.full_like( network.ext_t, network.params["stell_const_dc"]))
    network.ext_amp_intrnrn = h.Vector(np.full_like( network.ext_t, network.params["intrnrn_dc_amp"]))
    if network.params["input_id"]:
        add_extra_inputs(network)
    
    for stell in network.stellate_cells:
        instr_utils.set_stell_range_variables(stell,params)
        stell.ext_dc.dur = sim_dur
        network.ext_amp_stell.play(stell.ext_dc._ref_amp, network.ext_t, True)
        instr_utils.set_intial_noise(stell,params["stell_init_noise"],noise_seed=params["init_noise_seed"])
        # Noise
        instr_utils.set_noise(stell,params["stell_noise"],noise_seed=params["noise_seed"])
        #recorders
        instr_utils.setup_recorders(stell,params["record_handle_stell"],recorder_dt)


    for interneuron in network.interneurons:
        interneuron.ext_dc.dur = sim_dur
        network.ext_amp_intrnrn.play(interneuron.ext_dc._ref_amp, network.ext_t, True)
        instr_utils.set_intrnrn_range_variables(interneuron,params)
        instr_utils.set_intial_noise(interneuron,params["intrnrn_init_noise"],noise_seed=params["init_noise_seed"])
        instr_utils.set_noise(interneuron,params["intrnrn_noise"],noise_seed=params["noise_seed"])        
        #recorders
        instr_utils.setup_recorders(interneuron,params["record_handle_intrnrn"],recorder_dt)




def add_extra_inputs(network):
    if network.params["input_id"]=="pulse":
        pulse_width = network.params["extra_params"]["pulse_width"]
        amp = network.params["extra_params"]["amp"]
        start_0 = network.params["extra_params"]["start"]
        ipi = network.params["extra_params"]["ipi"]
        first_input = network.params["extra_params"]["first_cell_input"]
        first_input_gid = network.params["extra_params"]["first_cell_input"]+network.params["N_stell"]
        second_input_gid = (first_input+1)%2 + network.params["N_stell"]
        start_1 = (start_0 + ipi/2)

        pulse_0=instr_utils.generate_pulse_train(duration_ms=network.params["sim_dur"], pulse_width_ms=pulse_width, 
                                     ipi_ms=ipi, start_delay_ms=start_0, amplitude=amp)
        pulse_1=instr_utils.generate_pulse_train(duration_ms=network.params["sim_dur"], pulse_width_ms=pulse_width, 
                                     ipi_ms=ipi, start_delay_ms=start_1, amplitude=amp)
        
        pc = h.ParallelContext()
        if pc.gid_exists(first_input_gid):
            pc.gid2cell(first_input_gid).switch_pulse = h.Vector(pulse_0)
            pc.gid2cell(first_input_gid).switch_pulse_iclamp= h.IClamp(pc.gid2cell(first_input_gid).soma(0.5))
            pc.gid2cell(first_input_gid).switch_pulse_iclamp.dur = network.params["sim_dur"]
            pc.gid2cell(first_input_gid).switch_pulse.play(
                pc.gid2cell(first_input_gid).switch_pulse_iclamp._ref_amp, network.ext_t, True)

        if pc.gid_exists(second_input_gid):
            pc.gid2cell(second_input_gid).switch_pulse = h.Vector(pulse_1)
            pc.gid2cell(second_input_gid).switch_pulse_iclamp= h.IClamp(pc.gid2cell(second_input_gid).soma(0.5))
            pc.gid2cell(second_input_gid).switch_pulse_iclamp.dur = network.params["sim_dur"]
            pc.gid2cell(second_input_gid).switch_pulse.play(
                pc.gid2cell(second_input_gid).switch_pulse_iclamp._ref_amp, network.ext_t, True)
