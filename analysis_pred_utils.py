"""
Analysis functions for predictive coding simulations
"""
import numpy as np
from scipy import signal,stats
from scipy.signal import butter, sosfilt
import analysis_utils as a_utils

def separate_fields(spks):
    """Separates spike times into different fields based on a calculated threshold.
    
    Args:
        spks (list of float): A list of spike times.
    Returns:
        list of list of float: A list where each element is a list of spike times 
        that are categorized as the same field.
    """

    threshold=np.max(np.diff(np.sort(np.diff(spks))))
    separated_fields = [[]]
    for i,spk in enumerate(spks[:-1]):
        if (spks[i+1]-spk)>threshold:
           separated_fields[-1].append(spk)
           separated_fields.append([])
        else:
            separated_fields[-1].append(spk)
    return separated_fields

def calculate_field_positions(separate_fields_x,true_pos):
    """Calculate the positions at each spike based on true (interneuron) positions.

    Args:
        separate_fields_x (list of lists): A list where each element is a list 
            of field indices.
        true_pos (numpy.ndarray): A numpy array containing the true positions.
    Returns:
        list of lists: A list where each element is a list of true positions 
        corresponding to the field indices.
    """

    field_pos = []
    for field in separate_fields_x:
        field_pos.append(list(true_pos[np.floor(np.array(field)).astype("int")]))
    return field_pos

def find_rates_of_fields(separate_fields_x):
    """Calculate the instantaneous rates of each field.
    
    Args:
        separate_fields_x (list): A list of fields.
    Returns:
        list: A list of lists, where each inner list contains the rates 
        calculated for the corresponding field.
    """

    rates = []
    for i,field in enumerate(separate_fields_x):
        rates.append(list(a_utils.instant_rates(field,100)))
    return rates


def butter_lowpass(cutoff, fs, order=5):
    """Designs a lowpass Butterworth filter.
    
    Args:
        cutoff (float): The cutoff frequency of the filter.
        fs (float): The sampling frequency of the signal.
        order (int, optional): The order of the filter. Default is 5.
    Returns:
        ndarray: designed filter.
    """
    
    return butter(order, cutoff, fs=fs, btype='low', analog=False,output="sos")

def butter_lowpass_filter(data, cutoff, fs, order=5):
    """Apply a Butterworth lowpass filter to the given data.

    Args:
        data (array-like): The input signal data to be filtered.
        cutoff (float): The cutoff frequency of the filter.
        fs (float): The sampling frequency of the input signal.
        order (int, optional): The order of the filter. Default is 5.
    Returns:
        array-like: The filtered signal.
    """

    sos = butter_lowpass(cutoff, fs, order=order)
    y = sosfilt(sos, data)
    return y

def instant_rate_low_pass(stell_spikes_l,sim_dur):
    """Apply a low-pass filter to the instantaneous firing rates.
    
    This function bins the spike times and then applies a Butterworth low-pass filter
    at the cutoff frequency.

    Args:
        stell_spikes_l (list): List of spike times for stellate cells.
        sim_dur (float): Duration of the simulation in milliseconds.
    Returns:
        numpy.ndarray: The filtered instantaneous firing rate.
    """

    # Filter requirements.
    order = 6
    fs = 1000       # sample rate, Hz
    cutoff = 35  # desired cutoff frequency of the filter, Hz
    # Get the filter coefficients so we can check its frequency response.
    binned_spikes = a_utils.bin_spike_ms(stell_spikes_l, sim_dur)
    filtered = butter_lowpass_filter(binned_spikes, cutoff, fs, order)
    return filtered
    

def decode_pos_by_intrnrn(intrnrn_spks_l,params):
    """Decode position from interneuron spikes.

    Args:
        intrnrn_spks_l (list of lists): List of interneuron spikes.
        params (dict): Dictionary containing simulation parameters
    Returns:
        numpy.ndarray: Decoded positions.
    """
    
    N_per_sheet= params['N_per_sheet']
    n_phases= params['n_phases']
    sim_dur= params['sim_dur']
    lamb= params['lambda0']
    cell_phases = (np.arange(0,N_per_sheet)*(2*np.pi/n_phases))%(2*np.pi)
    t_stell=instant_rate_low_pass(intrnrn_spks_l,sim_dur)
    decoded_pos=((lamb/(2*np.pi))*((np.angle(np.sum((t_stell*np.exp(1j*cell_phases[:,np.newaxis])),axis=0)))))%lamb
    return decoded_pos



def bin_pos(field,period=None,res=0.098125):
    """Bin position array deriver from spikes.

    Args:
        field (array-like): Array of spike positions.
        period (float, optional): The period over which to bin the positions. Defaults to 2 * np.pi.
        res (float, optional): The resolution of the bins. Defaults to 0.098125.
    Returns:
        numpy.ndarray: Binned positions
    """
    
    if period is None:
        period = 2*np.pi
    pos_arr = np.arange(0,period,res)
    binned_pos= np.zeros_like(pos_arr)
    for spk in field:
        idx = (np.abs(pos_arr - spk)).argmin()
        binned_pos[idx]=1
    return binned_pos

def convole_field_pos(field_pos_x,params,win_size=6):
    """Convolves the input field positions derived from spikes to generate a firing field.
    
    Args:
        field_pos_x (list): List of field positions to be convolved.
        params (dict): Dictionary containing simulation parameters
        win_size (int, optional): Size of the Gaussian window. Default is 6.
    Returns:
        np.ndarray: Array of convolved field positions, truncated to the minimum length.
    """

    convoled = []
    win = signal.windows.hann(20)
    min_ax = 9999
    for field in field_pos_x:
        binned_pos = bin_pos(field,period=params["lambda0"],res=params["lambda0"]/params["n_phases"])
        win = signal.windows.gaussian(len(binned_pos),win_size,sym=False)
        filtered = circular_convolve(binned_pos,win)/sum(win)
        ax= len(filtered)
        if ax<min_ax:
            min_ax=ax
        convoled.append(filtered)
    for i,field in enumerate(convoled):
        convoled[i] = field[:min_ax]

    return np.roll(np.array(convoled),0,axis=1)
def remove_uneven_fields(cell_spks_x,bounds):
    """Removes spikes outside the field based on a threshold

    Args:
        cell_spks_x (array-like): Array of cell spike positions.
        bounds (tuple): A tuple containing the lower and upper 
            bounds of the direction change in the network.
    Returns:
        list: spikes times after removing escape spikes
    """
    
    cell_spks_x = cell_spks_x[np.array(cell_spks_x)>bounds[0]]
    cell_spks_x = cell_spks_x[np.array(cell_spks_x)<bounds[1]]
    separate_fields_x = separate_fields(cell_spks_x)
    lengths = []
    for field in separate_fields_x:
        lengths.append(len(field))
    threshold = np.max(np.diff(np.sort(lengths)))
    separated_fields = []
    for field in separate_fields_x:
        if len(field)>threshold:
            separated_fields.append(field)
    return separated_fields


def signed_arc_length(theta1, theta2, period=2*np.pi):
    """Calculate the signed arc length between two angles.
    
    Args:
        theta1 (float): The first angle in radians.
        theta2 (float): The second angle in radians.
        period (float, optional): The period of the circle, default is 2*pi.
    Returns:
        float: The signed arc length between theta1 and theta2. Positive if the 
        clockwise distance is shorter, negative if the counterclockwise distance 
        is shorter.
    """
    r=1
    theta1 = theta1 % period
    theta2 = theta2 % period
    delta_theta_cw = (theta1 - theta2) % period
    delta_theta_ccw = (theta2 - theta1) % period
    length_cw = r * delta_theta_cw
    length_ccw = r * delta_theta_ccw
    
    if length_cw<length_ccw:
        return length_cw
    else:
        return -length_ccw
        

def calc_field_size(stell_spks_l,true_pos,params,bounds):
    """Calculate the average field size w.r.t to positions.
    
    Args:
        stell_spks_l (list of lists): List of spike times for stellate cells.
        true_pos (array-like): Array of true positions.
        params (dict): Dictionary containing simulation parameters.
        bounds (tuple): Tuple containing the bounds for the analysis.

    Returns:
        float: The average field size.
    """

    field_sizes=[]
    for cell_to_analyse in range(params["N_per_sheet"]):
        cell_spks_r = np.array(stell_spks_l[cell_to_analyse])
        cell_spks_l = np.array(stell_spks_l[cell_to_analyse+params["N_per_sheet"]])
        separate_fields_r = remove_uneven_fields(cell_spks_r,(params["allothetic_dur"],bounds[0]))
        separate_fields_l = remove_uneven_fields(cell_spks_l,(bounds[0]+params["allothetic_dur"],bounds[1]))
        field_pos_l = calculate_field_positions(separate_fields_l,true_pos)
        field_pos_r = calculate_field_positions(separate_fields_r,true_pos)
        for field in field_pos_r:
            field_sizes.append(np.abs(circular_difference(field[-1], field[0],params["lambda0"])))
        for field in field_pos_l:
            field_sizes.append(np.abs(circular_difference(field[-1], field[0],params["lambda0"])))
    return np.mean(field_sizes)



def circular_difference(angle1, angle2,period=None):
    """Unsigned arc length

    Args:
        angle1 (float): The first angle in radians.
        angle2 (float): The second angle in radians.
        period (float, optional): The period of the circular scale. 
            Defaults to 2*pi.
    Returns:
        float: The smallest difference between the two angles.
    """

    if period is None:
        period = 2*np.pi
    diff = (angle2 - angle1) % period
    if diff > period/2:
        diff -= period
    return diff

def calc_speed_of_network(stell_spks_l,params,win_size=100):
    """Calculate the speed of the network based through decoded position.
    
    This function calculates the speed of the network by decoding the position 
    of the network and calculating the slope of the unwrapped position. 
    It assumes the input DC to the network was in the form of predictive coding simulations
    Same functionality as in analysis_utils.calc_speed_of_network().

    Args:
        stell_spks_l (list): List of spike times for stellate cells.
        params (dict): Dictionary containing simulation parameters
        win_size (int, optional): Window size for calculating instantaneous rates. Defaults to 100.
    Returns:
        float: The slope representing the speed of the network.
    """

    N_per_sheet= params['N_per_sheet']
    n_phases= params['n_phases']
    sim_dur= params['sim_dur']
    lamb = params['lambda0']
    allothetic_dur = params["allothetic_dur"]
    cell_phases = (np.arange(0,N_per_sheet)*(2*np.pi/n_phases))%(2*np.pi)
    cell_phases = np.concatenate((cell_phases,cell_phases))
    #instantaneous rates
    t_stell=a_utils.instant_rate_all(stell_spks_l,sim_dur,win_size)
    decoded_pos=((lamb/(2*np.pi))*((np.angle(np.sum((t_stell*np.exp(1j*cell_phases[:,np.newaxis])),axis=0)))))%lamb
    decoded_pos_unwrapped=np.unwrap(decoded_pos,period=lamb)
    x= (np.arange(0,params['sim_dur'])[:])/1000
    slope=stats.linregress(x[:params["extra_params"]["dir_change_t"]],np.unwrap(decoded_pos_unwrapped[:params["extra_params"]["dir_change_t"]])).slope
    return slope

def circular_convolve(x, h):
    """Convolution for circular data.
    
    Args:
        x (array-like): Data Array
        h (dict): Kernel
    Returns:
        numpy.ndarray: Convolved array.
    """

    N = len(x)
    if len(h) != N:
        raise ValueError("Inputs must have the same length")
    
    # Compute the circular convolution using FFT
    X = np.fft.fft(x)
    H = np.fft.fft(h)
    Y = X * H
    y = np.fft.ifft(Y)
    
    # Roll the result to correct the shift
    y = np.roll(y, -N//2)
    return np.real(y)

def calc_bias_stell_intrnrn(stell_decoded,intrnrn_decoded,params,bounds):
    """Calculate the positional bias between stellate and and interneuron.

    Here Interneuron is assumed to be the true position of the animal.
    
    Args:
        stell_decoded (np.ndarray): Decoded signal from stellate cells.
        intrnrn_decoded (np.ndarray): Decoded signal from interneurons.
        params (dict): Dictionary containing simulation parameters.
        bounds (tuple): Tuple containing the bounds (start, end) for the left 
            and right simulations.
    Returns:
        float: The calculated bias between the stellate and interneuron decoded signals.
    """

    period = params["lambda0"]
    allothetic_dur = params["allothetic_dur"]
    
    stell_decoded_r = stell_decoded[allothetic_dur:bounds[0]]
    intrnrn_decoded_r = intrnrn_decoded[allothetic_dur:bounds[0]]
    stell_decoded_l = stell_decoded[bounds[0]+allothetic_dur:bounds[1]]
    intrnrn_decoded_l = intrnrn_decoded[bounds[0]+allothetic_dur:bounds[1]]
    delta_theta_cw_r = (stell_decoded_r - intrnrn_decoded_r) % period
    delta_theta_ccw_r = (intrnrn_decoded_r - stell_decoded_r) % period
    r_bias = np.zeros_like(delta_theta_cw_r)
    r_bias[delta_theta_cw_r<delta_theta_ccw_r]=delta_theta_cw_r[delta_theta_cw_r<delta_theta_ccw_r]
    r_bias[delta_theta_cw_r>=delta_theta_ccw_r]=-delta_theta_ccw_r[delta_theta_cw_r>=delta_theta_ccw_r]
    delta_theta_cw_l = (stell_decoded_l - intrnrn_decoded_l) % period
    delta_theta_ccw_l = (intrnrn_decoded_l - stell_decoded_l) % period
    l_bias = np.zeros_like(delta_theta_cw_l)
    l_bias[delta_theta_cw_l<delta_theta_ccw_l]=delta_theta_cw_l[delta_theta_cw_l<delta_theta_ccw_l]
    l_bias[delta_theta_cw_l>=delta_theta_ccw_l]=-delta_theta_ccw_l[delta_theta_cw_l>=delta_theta_ccw_l]
    l_bias=-l_bias
    bias = np.mean(np.concatenate((r_bias,l_bias)))
    return bias
    
    
    
def calc_predictive_code(stell_spikes_l,intrnrn_spks_l,params,sim_num):
    """Calculate the predictive code from spikes

    This first decodes the position of stellates and interneurons and uses
    calc_bias_stell_intrnrn() to calculate the predictive code.
    
    Args:
        stell_decoded (np.ndarray): Decoded signal from stellate cells.
        intrnrn_decoded (np.ndarray): Decoded signal from interneurons.
        params (dict): Dictionary containing simulation parameters.
        bounds (tuple): Tuple containing the bounds (start, end) for the left 
            and right simulations.
    Returns:
        float: The calculated bias between the stellate and interneuron decoded signals.
    """

    bounds = (params["extra_params"]["dir_change_t"],int(params["sim_dur"]))
    true_pos = decode_pos_by_intrnrn(intrnrn_spks_l,params,win_size=40)
    stell_decoded = a_utils.decode_pos(stell_spikes_l,params)
    bias=calc_bias_stell_intrnrn(stell_decoded,true_pos,params,bounds)
    return bias
    
def calc_inhib_g_at_first_and_last_spike(stell_spks,stell_syn_inhib_g):
    """Get value of inhibitory synaptic conductance at first and last spike of all fields

    This function is used for plots in Figure 5.
    
    Args:
        x (array-like): Data Array
        h (dict): Kernel
    Returns:
        numpy.ndarray: Convolved array.
    """
    separate_fields = a_utils.separate_fields(stell_spks)
    fspk_inhib_g = []
    lspk_inhib_g = []
    for cell,fields in separate_fields.items():
        for field in fields:
                fspk = int(field[0])
                lspk = int(field[-1])
                fspk_inhib_g.append(stell_syn_inhib_g[cell][fspk])
                lspk_inhib_g.append(stell_syn_inhib_g[cell][lspk])
    return np.array([np.mean(fspk_inhib_g),np.mean(lspk_inhib_g)])