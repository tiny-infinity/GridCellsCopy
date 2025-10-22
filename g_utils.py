from scipy.ndimage import gaussian_filter1d
import sim_utils as s_utils
import analysis_utils as a_utils
from network_configs.instrumentations.trajectory1D import Trajectory1D
from scipy import stats

import matplotlib.pyplot as plt
import numpy as np

"""
New utility functions for analysis

"""
def decode_func(sim_id, sim_num, n_trials=10, sim_dur = float(60000)):
    """
    Take in Sim ID, Number

    Return arrays of mean and standard deviation of decoded trajectories
    
    """

    params = s_utils.load_sim_params(sim_id=sim_id)["0"]
    traj = Trajectory1D(params, save_mem=False)
    dc_inmput = traj.intrnrn_dc
    decoded_positions = np.zeros((n_trials, int(params["sim_dur"] - traj.init_allothetic_dur)))

    file_path_stell = f"data/{sim_id}/stell_spks_{sim_id}.hdf5"

    with h5py.File(file_path_stell, "r") as f:
        print("Top-level groups:", list(f.keys()))
        for grp in f.keys():
            print(grp, "contains:", list(f[grp].keys()))

    for tr in range(n_trials):
        print("Trial Number ", tr)
        stell_spks, _ = s_utils.load_spikes(sim_id, sim_num)
        decoded_positions[tr, :] = a_utils.decode_pos(stell_spks, params, win_size=40, t_start=int(traj.init_allothetic_dur))
        sim_num += 1

    mean_decode_posns = stats.circmean(decoded_positions, axis=0)
    std_decode_posns = stats.circstd(decoded_positions, axis=0)

    return mean_decode_posns, std_decode_posns


def posn_vel_input(sim_id, sim_num, n_trials=10, sim_dur = float(60000)):
    """
    Takes in Simulation ID, returns position, velocity and time arrays 
    
    """
    params = s_utils.load_sim_params(sim_id=sim_id)["0"]
    traj = Trajectory1D(params, save_mem=False)
    t_ms = np.arange(traj.init_allothetic_dur,params["sim_dur"])
    t_ms_idx = (t_ms/0.025).astype('int')
    position_input = traj.pos_input[t_ms_idx]
    velocity_input = traj.vel_input[t_ms_idx]
    t_s= np.linspace(traj.init_allothetic_dur/1000,(params['sim_dur']/1000),int(params["sim_dur"]-traj.init_allothetic_dur))

    return position_input, velocity_input, t_s


def plot_mult_stddev(sim_id_list,sim_num=0,save=None,num_trials=10):
    """
    Plots the standard deviation time series for a given simulation ID.
    Input:
        sim_id_list (list of strings) : list of Simulation IDs for which deviation is to be plotted
    
    
    """

    _,_,t_s=posn_vel_input(sim_id_list[0],sim_num)
    plt.figure(figsize=(22,10))
    for sim_id in sim_id_list:
        _,std = decode_func(sim_id,sim_num)
        plt.plot(t_s,std,label=f"{sim_id}")
    plt.legend()
    plt.title(f"Standard Deviation of decoded trajectories across {num_trials} trials")
    plt.ylabel("Standard Deviation (in rads)")
    plt.xlabel("Time (in s)")
    plt.show()
    if save:
        plt.savefig(f"{save}.png")
    return None

def plot_mult_error(sim_id_list,sim_num=0,save=None,num_trials=10):
    """
    Plots the circular error for all Sim IDs in the input list
    """
    pos_input,_,t_s=posn_vel_input(sim_id_list[0],sim_num)
    plt.figure(figsize=(22,10))
    for sim_id in sim_id_list:
        mean,_ = decode_func(sim_id,sim_num)
        error = circular_error(pos_input[:-1],mean[1:])
        plt.plot(t_s[1:],error,label=f"{sim_id}")
    plt.title(f"Position decoding error")
    plt.legend()
    plt.ylabel("Error (in rads)")
    plt.xlabel("Time (in s)")
    plt.show()
    if save:
        plt.savefig(f"{save}.png")
    return None

def com_raster(spk_array,convolve=True,gauss_window=10):
    """
    Finds the center of mass of neural activity using binned spike data
    Input:
        spk_array(ndarray) : matrix whose rows are cell indices and columns are millisecond timesteps
        gaussian
        convolve (boolean) : True if time series is to be smoothed out using a Gaussian filter
        gauss_window (float) : Standard deviation of Gaussian window
    Return:
        Time series of center of mass
    """
    com_time_series=[]
    for j in range(len(spk_array[0])):
        num,den=0,0
        for i in range(len(spk_array)):

            if spk_array[i][j]==1:
                num+=i
                den+=1
        if den==0:
            com_time_series.append(0)
        else:
            com_time_series.append(num/den)

    if convolve:
        return gaussian_filter1d(com_time_series,sigma=gauss_window)
    else:
        return com_time_series

def com_deviation(sim_id,num_sims=10,sim_dur_ms=60000):
    """
    Time series of standard deviation of center of mass of neural activity across multiple trials

    """
    mult_path=[]
    for i in range(num_sims):
        stell_spks,intrnrn_spks = s_utils.load_spikes(sim_id=sim_id,sim_num=i)
        binned_stell = a_utils.bin_spike_ms(stell_spks,sim_dur=sim_dur_ms)
        com_arr = com_raster(binned_stell)
        mult_path.append(com_arr)

    mult_arr = np.array(mult_path)

    return np.std(mult_arr,axis=0)


def circular_error(var1,var2):

    """Calculate the circular error between two angles.

    Parameters:
        var1 (float or array-like): The first angle in radians.
        var2 (float or array-like): The second angle in radians.
    Returns:
        numpy.ndarray: The minimum absolute difference between the angles
    """
    circum = 2*np.pi
    return np.abs(np.min(np.vstack((np.abs(var1-var2),circum-np.abs(var1-var2))),axis=0))