from scipy.ndimage import gaussian_filter1d
import sim_utils as s_utils
import analysis_utils as a_utils
from network_configs.instrumentations.trajectory1D import Trajectory1D
from scipy import stats
import h5py
import matplotlib.pyplot as plt
import numpy as np

"""
New utility functions for analysis

"""
def decode_func(sim_id, sim_num, n_trials, sim_dur = float(60000)):
    """
    Take in Sim ID, Number

    Return arrays of mean and standard deviation of decoded trajectories
    
    """
    print("Number of trials:" , n_trials)
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


def posn_vel_input(sim_id, sim_num, n_trials, sim_dur = float(60000)):
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


def plot_mult_stddev(sim_id_list,sim_num=0,save=None,num_trials=10,allo_time=None):
    """
    Plots the standard deviation time series for a given simulation ID.
    Input:
        sim_id_list (list of strings) : list of Simulation IDs for which deviation is to be plotted
    
    
    """
    print("num_trials : ", num_trials)
    _,_,t_s=posn_vel_input(sim_id_list[0],sim_num,n_trials=num_trials)
    plt.figure(figsize=(22,10))
    for sim_id in sim_id_list:
        _,std = decode_func(sim_id,sim_num,n_trials=num_trials)
        plt.plot(t_s,std,label=f"{sim_id}")
    plt.legend()
    plt.title(f"Standard Deviation of decoded trajectories across {num_trials} trials")
    plt.ylabel("Standard Deviation (in rads)")
    plt.xlabel("Time (in s)")
    
    if allo_time:
        plt.axvline(x=allo_time,ls='--',lw=1,color='black')
    
    if save:
        plt.savefig(f"{save}.png",dpi=300)
    
    plt.show()
    return None

def plot_mult_error(sim_id_list,sim_num=0,save=None,num_trials=10):
    """
    Plots the circular error for all Sim IDs in the input list
    """
    pos_input,_,t_s=posn_vel_input(sim_id_list[0],sim_num,n_trials=num_trials)
    plt.figure(figsize=(22,10))
    for sim_id in sim_id_list:
        mean,_ = decode_func(sim_id,sim_num,n_trials=num_trials)
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


def total_error(sim_id,sim_num=0,num_trials=10):
    pos_input,_,t_s=posn_vel_input(sim_id,sim_num,n_trials=num_trials)
    mean, std = decode_func(sim_id=sim_id,sim_num=sim_num,n_trials=num_trials)

    error = np.sum(circular_error(pos_input[:-1],mean[1:]))

    return error



def com_raster(spk_array,convolve=False,gauss_window=20):
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
    mean_time_series=[]
    var_time_series=[]
    for j in range(len(spk_array[0])): #for each timestep
        num,den,num2=0,0,0
        for i in range(len(spk_array)): #for each cell

            if spk_array[i][j]==1:
                num+=i
                den+=1
                num2+=i**2
        if den==0:
            mean_time_series.append(0)
            var_time_series.append(0)
        else:
            mean_time_series.append(num/den)
            var_time_series.append((num2/den) - (num/den)**2)


    if convolve:
        return gaussian_filter1d(mean_time_series,sigma=gauss_window), gaussian_filter1d(np.sqrt(var_time_series),sigma=gauss_window)
    else:
        return mean_time_series, np.sqrt(var_time_series)


def plot_mult_com(sim_id_list,sim_dur=60000,sim_num=0,convolve=False):
    """
    Plots the time series of the Center of Mass of spike activity for multiple simulations
    """
    time_arr = np.linspace(0,sim_dur,sim_dur)
    plt.figure(figsize=(22,12))
    for sim_id in sim_id_list:
        stell_spks,_ = s_utils.load_spikes(sim_id,sim_num)
        binned_spks = a_utils.bin_spike_ms(stell_spks,sim_dur)
        mean,std=com_raster(binned_spks,convolve=convolve)
        plt.plot(time_arr,mean,label=f'{sim_id}',linewidth=0.5)
    plt.title(f"Center of Mass trajectories")
    plt.xlabel("Time (in ms)")
    plt.ylabel("Position of activity center")
    plt.legend()
    plt.show()

    for sim_id in sim_id_list:
        stell_spks,_ = s_utils.load_spikes(sim_id,sim_num)
        binned_spks = a_utils.bin_spike_ms(stell_spks,sim_dur)
        mean,std=com_raster(binned_spks,convolve=convolve)
        plt.plot(time_arr,std,label=f'{sim_id}',linewidth=0.5)
    plt.title(f"Variance in trajectories")
    plt.xlabel("Time (in ms)")
    plt.ylabel("Variance")
    plt.legend()
    plt.show()

    plot_power_spectrum(com_raster(binned_spks),xmax=250)
        

def plot_power_spectrum(time_series,title=None,xmax=None):
    """
    Plots the power spectrum of any given time series
    
    """
    fft_freq,fft_sig,pow_spec = a_utils.calc_fft(time_series,T=0.001)
    plt.figure(figsize=(18,12))
    plt.plot(fft_freq,pow_spec)
    plt.xlabel("Frequency (in Hz)")
    plt.ylabel("Intensity")
    if xmax:
        plt.xlim(0,xmax)
    if title:
        plt.title(f"{title}")
    else:
        plt.title("Power Spectrum")

    plt.show()

    


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

def bin_spikes(sim_id,sim_num,sim_dur_ms):
    stell_spks,_ = s_utils.load_spikes(sim_id,sim_num)
    return a_utils.bin_spike_ms(stell_spks,sim_dur_ms)

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


def isi_distribution(spike_train,sim_dur):
    isi_list = []
    for i in range(len(spike_train)-1):
        isi_list.append((spike_train[i+1]-spike_train[i])/1000)
    isi_array = np.array(isi_list)
    plt.figure(figsize=(18,18))
    plt.hist(isi_array,bins=100,density=True)
    plt.show()