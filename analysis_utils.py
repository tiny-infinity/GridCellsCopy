import numpy as np
from scipy import signal,ndimage,integrate
import h5py
import os
import subprocess
import sim_utils as s_utils

def periodic_activity_all(stell_spikes_l, sim_dur, window_t, stdev):
    '''compute rates along the cells for all time window'''
    binned_spikes = bin_spike_ms(stell_spikes_l, sim_dur)
    binned_spikes_lim = np.array(
        np.hsplit(binned_spikes, sim_dur//window_t)).sum(axis=-1)
    kernel = signal.windows.gaussian(binned_spikes_lim.shape[1], stdev)
    filtered = ndimage.convolve1d(
        binned_spikes_lim, kernel, mode="wrap", axis=1)/np.sum(kernel)[None]
    return filtered

def bin_spike_ms(stell_spks_l:list, sim_dur:float)->np.ndarray:
    """Bin spike times into a binary matrix with millisecond resolution.
    
    Parameters
    ----------
    stell_spikes_l : list of lists  
        List of list where each sub-list corresponds to spike times of a cell.
    
    sim_dur : float or int
        Duration of the simulation in milliseconds.
    Returns
    -------
    t_stell : ndarray
        A binary matrix of shape (number of cells, sim_dur) where each row corresponds to 
        a cell and each column corresponds to a millisecond. A value of 1 indicates a spike 
        at that millisecond, and 0 indicates no spike.
    """
    
    t_stell = np.zeros((len(stell_spks_l), int(sim_dur)))
    for cell, stell_spks in enumerate(stell_spks_l):
        x = np.floor(stell_spks).astype("int")
        t_stell[cell, x] = 1
    return t_stell

def grid_field_sizes_neurons(stell_spks,sim_dur,avg=True,win_size_t=1000,win_size_n=3):
    period = grid_scale_neurons(stell_spks,sim_dur,win_size_t,win_size_n)
    cell_axis = np.arange(0,len(stell_spks))
    periodic_activity=periodic_activity_all(stell_spks,sim_dur,win_size_t,win_size_n)
    auc = []
    for win,periodic in enumerate(periodic_activity):
        peaks=signal.find_peaks(periodic)[0]
        periodic = periodic/np.max(periodic)
        for peak in peaks:
            start = int(cell_axis[peak] - period / 2)
            end = int(cell_axis[peak] + period / 2)
            if start < 0 or end >= cell_axis[-1]:
                continue  # Skip if the peak is near the signal boundaries
            auc.append(integrate_array(periodic[start:end]))
    if avg:
        return np.nanmean(auc)
    else:
        return auc

def grid_scale_neurons(stell_spikes_l, sim_dur, win_size=1000, stdev=3,avg=True):
    period_all = periodic_activity_all(
        stell_spikes_l, sim_dur, win_size, stdev)
    x = []
    if avg:
        for win in period_all:
            x.append(np.nanmean(np.diff(signal.find_peaks(win)[0])))
        return np.nanmean(x)
    else:
        for win in period_all:
            x.extend(np.diff(signal.find_peaks(win)[0]))
        return x

def integrate_array(arr, dx=1, axis=-1):
    return integrate.simpson(arr, dx=dx, axis=axis)


def build_and_return_matrix(sim_id:str=None,specs_file:str=None)->np.ndarray:
    """Builds and returns an connectivity matrix for a given simulation ID.
    
    Used to analyze the connectivity of a simulation.
    
    Parameters:
    -----------
        sim_id : str 
            The ID of the simulation to load parameters for.
        specs_file : str
            The path to a specifications file (not implemented).
    
    Returns:
    --------
        np.ndarray: The adjacency matrix generated for the given simulation.

    """
    
    if sim_id:
        params = s_utils.load_sim_params(sim_id)
        subprocess.run([f"python", f"network_configs/connections/{params['conn_id']}_config.py", "-i",f"{sim_id}"])
        with h5py.File(f"cache/matrix_{params['conn_id']}_{params['sim_id']}_{params['sim_num']}.hdf5", "r") as file:
            adj_matrix = np.array(file["matrix"])
        os.remove(f"cache/matrix_{params['conn_id']}_{params['sim_id']}_{params['sim_num']}.hdf5")
        return adj_matrix
    elif specs_file:
        raise NotImplementedError
    
def instant_rate_all(stell_spikes_l:list, sim_dur:float, stdev:float)->np.ndarray:
    """Calculate the instantaneous firing rate for all cells using Gaussian kernel convolution.
    
    Parameters:
        stell_spikes_l : list of lists  
            List of list where each sub-list corresponds to spike times of a cell.
        sim_dur : float 
            The duration of the simulation in milliseconds.
        stdev : float
            The standard deviation of the Gaussian kernel used for convolution.
    
    Returns:
        instantaneous rates : numpy.ndarray
            A 2D array (cell X t_ms) where each row corresponds to the instantaneous firing rate of a cell.
    """
    
    binned_spikes = bin_spike_ms(stell_spikes_l, sim_dur)
    kernal = np.vstack(
        (
            signal.windows.gaussian(len(binned_spikes[0]), stdev),
            np.zeros(len(binned_spikes[0])),
        )
    )
  
    filtered = signal.fftconvolve(binned_spikes, kernal, mode="same")
    return filtered

def spks_to_rate_reshaped(spks_l:list,params:dict,win_size:float=200)->np.ndarray:
    """Convert spike times to firing rates and reshape based on cell position.
    
    Parameters:
    -----------
        spks_l : list 
            List of list of spike times for each neuron.
        params : dict 
            Parameter dictionary containing simulation parameters.
        win_size : float, optional
            Window size for calculating the instantaneous rate. Default is 200.
    Returns:
    --------
        Instantaneous firing rates : np.ndarray
            Reshaped matrix of instantaneous firing rates.
    """
    
    inst_rate=instant_rate_all(spks_l,params['sim_dur'],win_size)
    inst_rate_reshaped=inst_rate.reshape(params['N_per_axis'],params['N_per_axis'],inst_rate.shape[1])
    inst_rate_reshaped = np.flip(inst_rate_reshaped,axis=0)
    return inst_rate_reshaped

def calc_fft(x,T = 0.025 * 1e-3):
    """Calculate the Fast Fourier Transform for a given signal.

    Parameters:
    -----------
        x : np.ndarray
            The signal to calculate the FFT for.
        T : float, optional
            The time step of the signal. Default is 0.000025s default dt for neuron.
    Returns:
    --------
        f : np.ndarray
            The frequency array.
        y : np.ndarray
            The FFT of the signal.
        power : np.ndarray
            The power spectrum of the signal.
    """

    from scipy.fft import fft, fftfreq
    x = signal.detrend(x)
    N = len(x)
    f = fftfreq(N, T)[: N // 2]
    y = fft(x)[: N // 2]
    power = (np.abs(y) ** 2)[: N // 2]
    return f, y, power

def decode_pos(stell_spikes_l,params,t_start=None,t_end=None,win_size=100):
    """Decode position for neuronal activity.
    
    Parameters:
    -----------
        stell_spikes_l : list of list
            List of list of spike times for each neuron.
        params : dict
            Parameter dictionary of simulation parameters.
        t_start : int, optional
            Start time for decoding. Default is 0.
        t_end : int, optional
            End time for decoding. Default is None.
        win_size : int, optional
            Window size for calculating the instantaneous rate. Default is 150.
    
    Returns:
    --------
        decoded_pos : np.ndarray
            The decoded position from.
    """
    N_per_sheet= params['N_per_sheet']
    n_phases= params['n_phases']
    sim_dur= params['sim_dur']
    lamb= params['lambda0']
    cell_phases = (np.arange(0,N_per_sheet)*(2*np.pi/n_phases))%(2*np.pi)
    cell_phases = np.concatenate((cell_phases,cell_phases))
    t_stell=instant_rate_all(stell_spikes_l,sim_dur,win_size)[:,t_start:t_end]
    #decode position
    decoded_pos=((lamb/(2*np.pi))*((np.angle(np.sum((t_stell*np.exp(1j*cell_phases[:,np.newaxis])),axis=0)))))%lamb
    return decoded_pos