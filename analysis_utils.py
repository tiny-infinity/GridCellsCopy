"""
Analysis functions
"""
import numpy as np
from scipy import signal,ndimage,integrate,stats
import h5py
import os
import subprocess
import sim_utils as s_utils
import copy

def periodic_activity_all(stell_spikes_l, sim_dur, window_t, stdev):
    """Instantaneous firing rate for all cells using Gaussian kernel convolution.
    
    Args:
        stell_spikes_l (list of lists): Spike times for stellate cells.
        sim_dur (int): Total duration of the simulation in milliseconds.
        window_t (int): Time window for binning spikes in milliseconds.
        stdev (float): Standard deviation for the Gaussian kernel used in convolution.
    Returns:
        numpy.ndarray: Filtered spike activity array with periodic activity highlighted.
    """

    binned_spikes = bin_spike_ms(stell_spikes_l, sim_dur)
    binned_spikes_lim = np.array(
        np.hsplit(binned_spikes, sim_dur//window_t)).sum(axis=-1)
    kernel = signal.windows.gaussian(binned_spikes_lim.shape[1], stdev)
    filtered = ndimage.convolve1d(
        binned_spikes_lim, kernel, mode="wrap", axis=1)/np.sum(kernel)[None]
    return filtered

def bin_spike_ms(stell_spks_l:list, sim_dur:float)->np.ndarray:
    """Bin spike times into a binary matrix with millisecond resolution.
    
    Args:
        stell_spikes_l (list of lists):Spike times for each neuron.
        sim_dur (float or int): Duration of the simulation in milliseconds.

    Returns:
        ndarray: A binary matrix of shape (number of cells, sim_dur) 
        where each row corresponds to a cell and each column corresponds to 
        a millisecond. A value of 1 indicates a spike at that millisecond, 
        and 0 indicates no spike.
    """
    
    t_stell = np.zeros((len(stell_spks_l), int(sim_dur)))
    for cell, stell_spks in enumerate(stell_spks_l):
        x = np.floor(stell_spks).astype("int")
        t_stell[cell, x] = 1
    return t_stell

def grid_field_sizes_neurons(stell_spks,sim_dur,avg=True,win_size_t=1000,win_size_n=3):
    """Calculate the grid field sizes along the neuron axis.

    Args:
        stell_spks (list of lists): Spike times for each neuron.
        sim_dur (float): Duration of the simulation in milliseconds.
        avg (bool, optional): If True, return the average grid field size. 
            If False, return all grid field sizes. Defaults to True.
        win_size_t (int, optional): Time window size in milliseconds for analysis. 
            Defaults to 1000.
        win_size_n (int, optional): Number of windows for analysis. Defaults to 3.

    Returns:
        float or list: Average grid field size if avg is True, otherwise a list 
        of all grid field sizes.
    """

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
    """Calculate the grid scale along the neuron axis.

    Args:
        stell_spikes_l (list of lists): Spike times for stellate neurons.
        sim_dur (int): Duration of the simulation in milliseconds.
        win_size (int, optional): Size of the window for periodic activity 
            calculation in milliseconds. Defaults to 1000.
        stdev (int, optional): Standard deviation for smoothing the periodic 
            activity. Defaults to 3.
        avg (bool, optional): If True, returns the average difference between peaks. 
            If False, returns all differences. Defaults to True.

    Returns:
        float or list: If `avg` is True, returns the average difference between peaks as a float.
        If `avg` is False, returns a list of differences between peaks.
    """

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
    """Integrate an array using Simpson's rule.
    """
    return integrate.simpson(arr, dx=dx, axis=axis)

def calc_speed_of_network(stell_spks_l,params,win_size=100):
    """Calculate the speed of the network.

    This function calculates the speed of the network by decoding the position 
    of the network and calculating the slope of the unwrapped position. 
    It assumes the input DC to network was constant
    
    Args:
        stell_spks_l (list): List of spike timings for stellate cells.
        params (dict): Dictionary containing the simulation parameters.
        win_size (int, optional): Window size for instantaneous rate calculation. 
            Defaults to 100.
    Returns:
        float: The slope representing the speed of the network.
    """

    N_per_sheet= params['N_per_sheet']
    n_phases= params['n_phases']
    sim_dur= params['sim_dur']
    lamb = (2*np.pi)
    cell_phases = (np.arange(0,N_per_sheet)*(2*np.pi/n_phases))%(2*np.pi)
    cell_phases = np.concatenate((cell_phases,cell_phases))
    #instantaneous rates
    t_stell=instant_rate_all(stell_spks_l[:],sim_dur,win_size)
    decoded_pos=((lamb/(2*np.pi))*((np.angle(np.sum((t_stell*np.exp(1j*cell_phases[:,np.newaxis])),axis=0)))))%lamb
    decoded_pos_unwrapped=np.unwrap(decoded_pos,period=lamb)
    x= (np.arange(0,params['sim_dur'])[:])/1000
    slope=stats.linregress(x[:],np.unwrap(decoded_pos_unwrapped[:])).slope
    return slope


def calc_firing_rates(spk_list,sim_dur):
    """Calculate the average firing rates of neurons.

    Args:
        spk_list (list of lists): Spike times of each neuron
        sim_dur (float): The duration of the simulation in milliseconds.
    Returns:
        float: The average firing rate of the neurons in Hz.
    """

    avg_fr_rates = []
    for cell in spk_list:
        avg_fr_rates.append(len(cell)/(sim_dur/1000))
    return np.mean(avg_fr_rates) #in Hz

def build_and_return_matrix(sim_id:str=None,specs_file:str=None)->np.ndarray:
    """Builds and returns an connectivity matrix for a given simulation ID.
    
    Used to analyze the connectivity of a simulation.
    
    Args:
        sim_id (str): The ID of the simulation to load parameters for.
        specs_file (str): The path to a specifications file (not implemented).
    
    Returns:
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
    
    Args:
        stell_spikes_l (list of lists): Spike times for each neuron.
        sim_dur (float): The duration of the simulation in milliseconds.
        stdev (float): The standard deviation of the Gaussian kernel used for convolution.
    
    Returns:
        numpy.ndarray: A 2D array (cell X t_ms) where each 
        row corresponds to the instantaneous firing rate of a cell.
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

def instant_rate(spike_train, sim_dur, stdev):
    """Calculate the instantaneous firing rate from a spike train using Gaussian convolution.

    Args:
        spike_train (list or array-like): A list or array of spike times.
        sim_dur (float): The duration of the simulation in milliseconds.
        stdev (float): The standard deviation of the Gaussian window.
    Returns:
        numpy.ndarray: The instantaneous firing rate as a numpy array.
    """

    t_stell = np.zeros(int(sim_dur))
    for i in spike_train:
        t_stell[int(np.floor(i))] = 1
    win = signal.windows.gaussian(len(t_stell), stdev)

    filtered = signal.convolve(
        t_stell, win, mode="same", method="fft") / sum(win)
    return filtered

def grid_props_2D(inst_rates_reshaped,t=-1000):
    """Calculate 2D grid properties.
    
    Calculate the grid score, grid scale, and grid size from the 2D autocorrelation of firing rate map.

    Args:
        inst_rates_reshaped (numpy.ndarray): 3D array of reshaped instantaneous 
            rates (cellxcellxtime).
        t (int): Time index to use for the calculation. Default is -1000.
    Returns:
        tuple: A tuple containing:
            - ``auto_corr`` (numpy.ndarray): Normalized 2D autocorrelation matrix.
            - ``grid_score`` (float): Grid score calculated from the autocorrelation matrix.
            - ``grid_scale`` (float): Median distance of the closest objects from the center.
            - ``grid_size`` (float): Median size of the grid cells.
    """

    auto_corr=signal.correlate2d(inst_rates_reshaped[:,:,t],inst_rates_reshaped[:,:,t],mode='full',boundary='wrap')
    auto_corr_mod = auto_corr.copy()
    auto_corr_mod[auto_corr_mod<(0.25*np.max(auto_corr))]=0
    features,nclusters=ndimage.label(auto_corr_mod,np.array([[1,1,1],[1,1,1],[1,1,1]]))
    obj = ndimage.find_objects(features)
    obj_dist_from_centre = np.zeros(nclusters)
    grid_size=[]
    center_object_idx = np.array([len(features)//2,len(features)//2])
    for i in range(nclusters):
        xdist=obj[i][0].stop-obj[i][0].start
        ydist=obj[i][1].stop-obj[i][1].start
        grid_size.append(np.mean([xdist,ydist]))
        row_cent_idx=int((obj[i][0].stop+obj[i][0].start)/2)
        col_cent_idx=int((obj[i][1].stop+obj[i][1].start)/2)
        obj_idx = np.array([row_cent_idx,col_cent_idx])
        obj_dist_from_centre[i]=np.linalg.norm(np.abs(obj_idx-center_object_idx))


    closest_obj_idx = np.argsort(obj_dist_from_centre)[1:7] #obj idx
    closest_obj = closest_obj_idx+1 #object number
    closest_obj_lims = np.zeros((len(closest_obj_idx),2,2))
    for i,idx in enumerate(closest_obj_idx):
        closest_obj_lims[i,:,0]=obj[idx][0].start,obj[idx][0].stop
        closest_obj_lims[i,:,1]=obj[idx][1].start,obj[idx][1].stop
    try:
        xmin,xmax=int(np.min(closest_obj_lims[:,:,0])),int(np.max(closest_obj_lims[:,:,0]))
        ymin,ymax=int(np.min(closest_obj_lims[:,:,1])),int(np.max(closest_obj_lims[:,:,1]))
    
        auto_corr_cropped=auto_corr_mod[xmin:xmax,ymin:ymax]
        rotations = [30,60,90,120,150]
        auto_corr_cropped_rotated = {}
        for i in rotations:
            auto_corr_cropped_rotated[i]=ndimage.rotate(auto_corr_cropped,i,reshape = False)

        correlations = {}
        for key,val in auto_corr_cropped_rotated.items():  
            correlations[key]=(stats.pearsonr(auto_corr_cropped.flatten(),val.flatten()))[0]

        min_corr = min(correlations[60],correlations[120])
        max_corr = max(correlations[30],correlations[90],correlations[150])
        grid_score= min_corr-max_corr
        grid_scale = np.median(np.sort(obj_dist_from_centre)[1:7])
        grid_size = np.median(grid_size)
        return  auto_corr/np.max(auto_corr),np.round(grid_score,2),round(grid_scale,2),round(grid_size,2)
    except:
        return auto_corr/np.max(auto_corr),np.nan,np.nan,np.nan

def spks_to_rate_reshaped(spks_l:list,params:dict,win_size:float=200)->np.ndarray:
    """Convert spike times to firing rates and reshape based on cell position.
    
    Args:
        spks_l (list of lists): Spike times for each neuron.
        params (dict): Parameter dictionary containing simulation parameters.
        win_size (float, optional): Window size for calculating the 
            instantaneous rate. Default is 200.
    Returns:
        np.ndarray: Reshaped matrix of instantaneous firing rates.
    """
    
    inst_rate = instant_rate_all(spks_l, params['sim_dur'], win_size)
    inst_rate_reshaped=inst_rate.reshape(params['N_per_axis'],params['N_per_axis'],inst_rate.shape[1])
    inst_rate_reshaped = np.flip(inst_rate_reshaped,axis=0)
    return inst_rate_reshaped

def calc_fft(x,T = 0.025 * 1e-3):
    """Calculate the Fast Fourier Transform for a given signal.

    Args:
        x (np.ndarray): The signal to calculate the FFT for.
        T (float, optional): The time step of the signal. Default is 0.000025s default dt for neuron.
    Returns:
        tuple: A tuple containing:
            - f (np.ndarray): The frequency array.
            - y (np.ndarray): The FFT of the signal.
            - power (np.ndarray): The power spectrum of the signal.
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
    
    Args:
        stell_spikes_l (list of lists): Spike times for each neuron.
        params (dict): Parameter dictionary of simulation parameters.
        t_start (int, optional): Start time for decoding. Default is 0.
        t_end (int, optional): End time for decoding. Default is None.
        win_size (int, optional): Window size for calculating the instantaneous 
            rate. Default is 150.
    Returns:
        np.ndarray: The decoded position
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

def clean_spikes(stell_spikes_l,order=1):
    """Cleans spikes from the given spike trains using a threshold-based approach.

    Args:
        stell_spikes_l (list of lists): Spike times for each neuron.
        order (float, optional): Factor to determine the strictness of threshold.
            
    Returns:
        list: A list of spike trains with cleaned spikes, where each spike train 
        is represented as a list of spike times.

    """
    stell_spks_clean = copy.deepcopy(stell_spikes_l)
    thresholds = []
    for i, stell in enumerate(stell_spikes_l):
        if len(stell) > 5:
            prev_spike = stell[0]
            sorted_isi = np.sort(np.diff(stell))
            thresholds.append(
                np.mean(sorted_isi[np.argmax(
                    np.diff(np.diff(sorted_isi))) + 2:])
            )  # from largest change in slope of ISI to the largest ISI
    threshold = order*np.abs(np.mean(thresholds))

    for i, stell in enumerate(stell_spikes_l):
        if len(stell) > 5:
            for j in range(len(stell) - 2):
                spk = stell[j + 1]
                prev_spike = stell[j]
                next_spike = stell[j + 2]
                delta_minus = spk - prev_spike
                delta_plus = next_spike - spk
                if (delta_minus + delta_plus > (threshold)) and abs(
                    delta_plus - delta_minus
                ) < threshold: #remove spikes that are in the middle of two fields
                    stell_spks_clean[i][j] = 0
    res_spks_clean = [[] for x in range(len(stell_spikes_l))]

    for i, stell in enumerate(stell_spks_clean):
        for spks in stell:
            if spks != 0:
                res_spks_clean[i].append(spks)
    return res_spks_clean



def separate_fields(stell_spikes_l,order=1):
    """
    Separates spike trains into grid fields based on a threshold.

    Args:
        stell_spikes_l (list of lists): Spike times for each neuron.

    Returns:
        dict: A dictionary with keys as cell indices and values as 
        lists of lists of separated grid fields,
    """
    grid_fields_all = {i: None for i in range(len(stell_spikes_l))}
    stell_spikes_l = clean_spikes(stell_spikes_l)
    for i, stell in enumerate(stell_spikes_l):
        if len(stell) > 5:
            prev_spike = stell[0]
            fields = [[prev_spike]]
            sorted_isi = np.sort(np.diff(stell))
            idx = np.argmax(np.diff(sorted_isi)) + 1
            threshold = (sorted_isi[idx] - (sorted_isi[idx] % 10))*order
            for j, spk in enumerate(stell[1:]):
                delta = spk - prev_spike
                if abs(delta) < abs(threshold):
                    fields[-1].append(spk)
                    prev_spike = spk
                else:
                    fields.append([spk])
                    prev_spike = spk
            grid_fields_all[i] = fields
    return grid_fields_all

def calc_grid_field_sizes_time(stell_spikes_l, avg=True):
    """Calculate the sizes of grid fields along the time axis.

    Args:
        stell_spikes_l (list of lists): Spike times for each neuron.
        avg (bool, optional): Whether to return the average field size. Defaults to True.

    Returns:
        float or list: If avg is True, returns the median field size. 
        Otherwise, returns a list of all field sizes.
    """
    grid_fields_all = separate_fields(stell_spikes_l)
    all_field_size = []
    for cell, fields in grid_fields_all.items():
        if fields != None:
            cell_field_size = []
            for a_field in fields:
                if len(a_field) > 3:
                    cell_field_size.append(a_field[-1] - a_field[0])
            all_field_size.extend(cell_field_size[1:-1])
    if avg:
        return np.median(all_field_size)
    else:
        return all_field_size

def calc_grid_scales_time(stell_spikes,avg=True):
    """Calculate the scale of grid fields along the time axis.

    Args:
        stell_spikes_l (list of lists): Spike times for each neuron.
        avg (bool, optional): Whether to return the average scale size. Defaults to True.

    Returns:
        float or list: If avg is True, returns the median scale size. Otherwise, 
        returns a list of all scale sizes.
    """    
    grid_fields_all = separate_fields(stell_spikes)
    all_cell_scales = []
    for cell, fields in grid_fields_all.items():
        if fields != None:
            cell_field_scales = []
            for a_field in fields:
                if len(a_field) > 3:
                    cell_field_scales.append(np.mean(a_field))

            all_cell_scales.extend(np.diff(cell_field_scales)[1:-1])

    # returns nan if number of fields too low
    if avg:
        return np.median(all_cell_scales)
    else:
        return all_cell_scales

def shift_fields_to_center(stell_spikes):
    """Shift the separated grid fields of spike trains to center them around zero.

    Args:
        stell_spikes (list of lists): Spike times for each neuron.

    Returns:
        dict: A dictionary with keys as cell indices and values as lists of 
        shifted grid fields.
    """

    separated_fields = separate_fields(stell_spikes)
    shifted_fields = {}.fromkeys(separated_fields.keys())
    for cell, fields in separated_fields.items():
        if fields != None:
            shifted_field_cell = []
            for a_field in fields[1:-1]:
                if len(a_field) > 0:
                    field_center = np.median(a_field)
                    shifted_field_cell.append(np.array(a_field) - field_center)
            shifted_fields[cell] = shifted_field_cell
    return shifted_fields



def generate_2d_video(sim_id,sheet_to_save=0):
    """Generates a video for the 2D model

    Video is saved as data/{sim_id}/{sim_id}_{sheet_to_save}.mp4

    Args:
        sim_id (str): The simulation ID to load spikes and parameters.
        sheet_to_save (int, optional): The sheet index to save (1-4). Defaults to 0. 
            Index 4 is the interneuron sheet.
    Returns:
        None
    Example:
        generate_2d_video("BaseModel2D", sheet_to_save=2)
    """
    
    import matplotlib.pyplot as plt
    import time
    import matplotlib.animation as animation
    import seaborn as sns
    plt.style.use('analysis/config/paper.mplstyle')
    plt.rcParams["ytick.labelsize"] = 22
    plt.rcParams["xtick.labelsize"] = 22
    
    t1 = time.perf_counter()
    stell_spikes_l,intrnrn_spikes_l=s_utils.load_spikes(sim_id=sim_id)
    params=s_utils.load_sim_params(sim_id=sim_id)   
    if list(params.keys())[0] == '0':
        params = params['0']
    n_per_sheet=params["N_per_sheet"]
    idx = np.full((4,2),n_per_sheet)*np.array([[0,1],[1,2],[2,3],[3,4]])
    print("Reshaping array")
    if sheet_to_save <=3:
        stell_spikes_arr =spks_to_rate_reshaped(
            stell_spikes_l[idx[sheet_to_save][0]:idx[sheet_to_save][1]], params, win_size=100)
    else:
        stell_spikes_arr = spks_to_rate_reshaped(
            intrnrn_spikes_l, params, win_size=100)
    stell_spikes_arr= np.flip(stell_spikes_arr,axis=0) # let matplotlib handle origin
    print(f"Generating video for {sim_id}, sheet: {sheet_to_save}")
    fig,axs= plt.subplots(1,1)
    fig.set_tight_layout(True)
    plot = plt.imshow(stell_spikes_arr[:,:,0],cmap=sns.cubehelix_palette(as_cmap=True,reverse=True),origin="lower",
                      vmax=np.max(stell_spikes_arr[:,:,int(params["stell_init_noise"][0]+500):int(params["sim_dur"])]))
    
    axs.set(xticks=np.arange(1,params["N_per_axis"],20),xticklabels=np.arange(1,params["N_per_axis"],20),
            yticks=np.arange(1,params["N_per_axis"],20),yticklabels=np.arange(1,params["N_per_axis"],20))
    axs.tick_params(which="both",direction="out",pad=0.5)
    axs.grid(False)
    def animate(i):
        stell_spikes_i = stell_spikes_arr[:,:,i]
        plot.set_data(stell_spikes_i)
        return [plot,]

    fps = 30 #higher fps higher temporal resolution (longer encoding time, displays might not support)
    total_samples = int(fps*(params['sim_dur']/1000))
    plot_frames =np.linspace(0,params['sim_dur']-1,total_samples,dtype='int')
    anim = animation.FuncAnimation(fig, animate,frames = plot_frames, interval =1, blit = True)
    writer = animation.FFMpegWriter(fps=fps,bitrate=8000) #bitrate: quality vs file_size
    
    anim.save(f"data/{params['sim_id']}/{params['sim_id']}_{sheet_to_save}.mp4", writer=writer)
    plt.close(fig)
    print(f"Video saved in data/{params['sim_id']}/{params['sim_id']}_{sheet_to_save}.mp4 ",
          "t_total (s):",round(time.perf_counter()-t1,2))

   