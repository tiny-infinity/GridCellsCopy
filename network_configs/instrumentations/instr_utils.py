"""
Utility functions for setting up recorders, noise, global variables and range 
variables in the network.
"""
from neuron import h
import numpy as np

def set_intial_noise(cell,noise_params,noise_seed=40):

    noise_seed = cell._gid+noise_seed
    cell.init_noise.dur = noise_params[0]
    cell.init_noise_t = np.arange(0, cell.init_noise.dur, h.dt)
    cell.init_noise_t_amp = h.Vector(
        np.random.default_rng(noise_seed).normal(
            noise_params[1],
            noise_params[2],
            cell.init_noise_t.shape,
        )
    )
    cell.init_noise_t = h.Vector(cell.init_noise_t)
    cell.init_noise_t_amp.play(
        cell.init_noise._ref_amp, cell.init_noise_t, True
    )

def set_noise(cell,noise_params,noise_seed=80):
    noise_seed = cell._gid+noise_seed
    cell.noise.dur = noise_params[0]
    cell.noise_t = np.arange(
        cell.init_noise.dur, cell.noise.dur, h.dt)
    cell.noise_amp = h.Vector(
        np.random.default_rng(noise_seed).normal(
            noise_params[1],
            noise_params[2],
            cell.noise_t.shape,
        )
    )
    cell.noise_t = h.Vector(cell.noise_t)
    cell.noise_amp.play(cell.noise._ref_amp, cell.noise_t, True)

def recursive_getattr(obj, attr_string):
    """Recursively retrieves an attribute or calls a method on an object based on a colon-separated string.

    Used to add recorders.

    Parameters:
        obj
            The object from which to retrieve the attribute or method.
        attr_string : str
            A colon-separated string representing the attribute or method to retrieve or call.
            If the string ends with ')', it indicates a method call.

    Returns:
        object: The final attribute or the result of the method call.
    """
    parts = attr_string.split(':')
    for part in parts:
        if part.endswith(")"):
            func_name = part[:part.find("(")]
            func = getattr(obj, func_name)
            obj = func()
        else:
            obj = getattr(obj, part)
    return obj
def setup_recorders(cell,recorder_handle,recorder_dt):
    """Sets up recorders for a given cell based on the provided recorder handle.

    For each parameter, create a vector to record from the given location in the cell obj.
    For spikes, create a spike detector and record the spikes.

    Parameters:
        cell
            The cell object for which the recorders are being set up.
        recorder_handle : dict 
            A dictionary containing recorder configurations.
        recorder_dt : float
            The time interval for recording data.
    Returns:
        object: The cell object with the recorders set up.

    """
    #dictionary poiting to the location of the parameter in the cell
    for record in recorder_handle:
        #Dealing with non-spike data
        if not record.endswith("spks") and recorder_handle[record]['state']==True and \
        (recorder_handle[record]['cells_to_record']=='all' or \
            cell._gid in list(recorder_handle[f'{record}']['cells_to_record'])):
            cell.recorder[f'{record}'] = h.Vector().record( \
                recursive_getattr(cell,recorder_handle[record]["loc"]),recorder_dt)
        
        #Dealing with spikes
        elif record.endswith("spks") and recorder_handle[record]['state']==True and \
            (recorder_handle[record]['cells_to_record']=='all' or \
             cell._gid in list(recorder_handle[f'{record}']['cells_to_record'])):
            
            cell.recorder[f'{record}'] = h.Vector()
            cell._spike_detector.record(cell.recorder[f'{record}'])
            
    return cell

def set_global_variables(params):
    #global variables
    h.hf_tau_input_stellate_mech = params['hf_tau']
    h.hs_tau_input_stellate_mech = params['hs_tau']
    h.phi_i_theta = params['phi_i_theta']
    h.omega_i_theta = params['omega_i_theta']
    
def set_stell_range_variables(stell,params):
    stell.soma(0.5).stellate_mech.ghbar = params['g_h_bar']
    stell.soma(0.5).stellate_mech.gnap_bar = params['gnap_bar']
    stell.soma(0.5).i_theta_stell.Amp = params['stell_theta_Amp']
    stell.soma(0.5).i_theta_stell.omega = params['stell_theta_omega']

def set_intrnrn_range_variables(interneuron,params):
    interneuron.soma(0.5).i_theta.Amp =params["Amp_i_theta"]
    