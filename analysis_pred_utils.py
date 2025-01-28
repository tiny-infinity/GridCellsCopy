"""
Analysis functions for predictive coding
"""
import numpy as np
from scipy import signal,stats
from scipy.signal import butter, sosfilt
import analysis_utils as a_utils

def separate_fields(spks):
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
    field_pos = []
    for field in separate_fields_x:
        field_pos.append(list(true_pos[np.floor(np.array(field)).astype("int")]))
    return field_pos

def find_rates_of_fields(separate_fields_x):
    rates = []
    for i,field in enumerate(separate_fields_x):
        rates.append(list(a_utils.instant_rates(field,100)))
    return rates


def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype='low', analog=False,output="sos")

def butter_lowpass_filter(data, cutoff, fs, order=5):
    sos = butter_lowpass(cutoff, fs, order=order)
    y = sosfilt(sos, data)
    return y

def instant_rate_low_pass(stell_spikes_l,sim_dur):
    # Filter requirements.
    order = 6
    fs = 1000       # sample rate, Hz
    cutoff = 35  # desired cutoff frequency of the filter, Hz
    # Get the filter coefficients so we can check its frequency response.
    binned_spikes = a_utils.bin_spike_ms(stell_spikes_l, sim_dur)
    filtered = butter_lowpass_filter(binned_spikes, cutoff, fs, order)
    

    return filtered
    

def decode_pos_by_intrnrn(intrnrn_spks_l,params,win_size=150):
    N_per_sheet= params['N_per_sheet']
    n_phases= params['n_phases']
    sim_dur= params['sim_dur']
    lamb= params['lambda0']
    cell_phases = (np.arange(0,N_per_sheet)*(2*np.pi/n_phases))%(2*np.pi)
    t_stell=instant_rate_low_pass(intrnrn_spks_l,sim_dur)
    decoded_pos=((lamb/(2*np.pi))*((np.angle(np.sum((t_stell*np.exp(1j*cell_phases[:,np.newaxis])),axis=0)))))%lamb
    return decoded_pos



def bin_pos(field,period=None,res=0.098125):
    if period is None:
        period = 2*np.pi
    pos_arr = np.arange(0,period,res)
    binned_pos= np.zeros_like(pos_arr)
    for spk in field:
        idx = (np.abs(pos_arr - spk)).argmin()
        binned_pos[idx]=1
    return binned_pos

def convole_field_pos(field_pos_x,params,win_size=6):
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
def signed_arc_length(theta1, theta2,period=2*np.pi):
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
def calc_bias_all_cells(stell_spks_l,true_pos,params,bounds,win_size=6,allothetic_dur=3000,ret_direction=False):
    sim_dur = params["sim_dur"]
    bias = []
    direction = []


    for cell_to_analyse in range(params["N_per_sheet"]):
        cell_spks_r = np.array(stell_spks_l[cell_to_analyse])
        cell_spks_l = np.array(stell_spks_l[cell_to_analyse+params["N_per_sheet"]])
        separate_fields_r = remove_uneven_fields(cell_spks_r,(allothetic_dur,bounds[0]-1000))
        separate_fields_l = remove_uneven_fields(cell_spks_l,(bounds[0]+allothetic_dur,bounds[1]))
        field_pos_l = calculate_field_positions(separate_fields_l,true_pos)
        field_pos_r = calculate_field_positions(separate_fields_r,true_pos)
        convoled_l = convole_field_pos(field_pos_l,params,win_size=win_size)
        convoled_l=np.mean(convoled_l,axis=0)
        field_pos_r = calculate_field_positions(separate_fields_r,true_pos)
        convoled_r = convole_field_pos(field_pos_r,params,win_size=win_size)
        convoled_r=np.mean(convoled_r,axis=0)
        pos = np.arange(0,2*np.pi,0.1)
        pos = np.arange(0,params["lambda0"],params["lambda0"]/params["n_phases"])
        try:

            means = np.mean(np.row_stack((convoled_l,convoled_r)),axis=0)
            peak_l=pos[np.argmax(convoled_l)]
            peak_r=pos[np.argmax(convoled_r)]
            dist=signed_arc_length(peak_l,peak_r,period = params["lambda0"])
            if dist>0:
                direction.append(1)
            else:
                direction.append(-1)

            peak_mean=pos[np.argmax(means)]
        except:
            continue
        
            # bias.append(np.abs(peak_mean-peak_r))
            # bias.append(np.abs(peak_mean-peak_l))
        bias.append(np.abs(circular_difference(peak_mean, peak_l,params["lambda0"])))
        bias.append(np.abs(circular_difference(peak_mean, peak_r,params["lambda0"])))

            # bias.append(np.abs(circular_difference(peak_l, peak_r)))

            
    # plt.figure()
    # plt.hist(bias)
    if not ret_direction:
        return np.mean(bias)
    else:
        if np.mean(direction)>0:
            return np.mean(bias),1
        else:
            return np.mean(bias),-1
def calc_bias_all_cells_first_spk(stell_spks_l,true_pos,params,bounds,win_size=6):
    #w.r.t to start of the field
    sim_dur = params["sim_dur"]
    allothetic_dur= params["allothetic_dur"]
    bias = []
    for cell_to_analyse in range(params["N_per_sheet"]):
        cell_spks_r = np.array(stell_spks_l[cell_to_analyse])
        cell_spks_l = np.array(stell_spks_l[cell_to_analyse+params["N_per_sheet"]])
        separate_fields_r = remove_uneven_fields(cell_spks_r,(allothetic_dur,bounds[0]))
        separate_fields_l = remove_uneven_fields(cell_spks_l,(bounds[0]+allothetic_dur,bounds[1]))
        field_pos_l = calculate_field_positions(separate_fields_l,true_pos)
        field_pos_r = calculate_field_positions(separate_fields_r,true_pos)
        convoled_l = convole_field_pos(field_pos_l,params,win_size=win_size)
        convoled_l=np.mean(convoled_l,axis=0)
        convoled_r = convole_field_pos(field_pos_r,params,win_size=win_size)
        convoled_r=np.mean(convoled_r,axis=0)
        pos = np.arange(0,params["lambda0"],params["lambda0"]/params["n_phases"])
        try:
            means = np.mean(np.row_stack((convoled_l,convoled_r)),axis=0)
            # plt.figure()
            # plt.plot(means)
            # plt.show()
            # peak_mean=pos[signal.find_peaks(means)[0][0]]
            peak_mean=pos[np.argmax(means)]

        except:
            continue
        for field in field_pos_r:
            bias.append(np.abs(circular_difference(peak_mean, field[0],params["lambda0"])))
        for field in field_pos_l:
            bias.append(np.abs(circular_difference(peak_mean, field[-1],params["lambda0"])))
    return np.mean(bias)
def calc_field_size(stell_spks_l,true_pos,params,bounds):

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
    """
    Calculate the circular difference between two angles in degrees.
    """
    if period is None:
        period = 2*np.pi
    diff = (angle2 - angle1) % period
    if diff > period/2:
        diff -= period
    return diff

def calc_speed_of_network(stell_spks_l,params,win_size=100):
    N_per_sheet= params['N_per_sheet']
    n_phases= params['n_phases']
    sim_dur= params['sim_dur']
    lamb = (2*np.pi) #params['n_phases']
    allothetic_dur = 3000 #when does the network stabilize?
    cell_phases = (np.arange(0,N_per_sheet)*(2*np.pi/n_phases))%(2*np.pi)
    # print(cell_phases)
    cell_phases = np.concatenate((cell_phases,cell_phases))
    #instantaneous rates
    t_stell=a_utils.instant_rate_all(stell_spks_l,sim_dur,win_size)
    decoded_pos=((lamb/(2*np.pi))*((np.angle(np.sum((t_stell*np.exp(1j*cell_phases[:,np.newaxis])),axis=0)))))%lamb
    decoded_pos_unwrapped=np.unwrap(decoded_pos,period=lamb)
    x= (np.arange(0,params['sim_dur'])[:])/1000
    slope=stats.linregress(x[:params["extra_params"]["dir_change_t"]],np.unwrap(decoded_pos_unwrapped[:params["extra_params"]["dir_change_t"]])).slope
    return slope

def calc_true_pos(stell_spks_l,params,win_size=100):
    #calculate speed of network
    N_per_sheet= params['N_per_sheet']
    n_phases= params['n_phases']
    sim_dur= params['sim_dur']
    lamb = (2*np.pi) #params['n_phases']
    allothetic_dur = 3000 #when does the network stabilize?
    cell_phases = (np.arange(0,N_per_sheet)*(2*np.pi/n_phases))%(2*np.pi)
    # print(cell_phases)
    cell_phases = np.concatenate((cell_phases,cell_phases))
    #instantaneous rates
    t_stell=a_utils.instant_rate_all(stell_spks_l,sim_dur,win_size)
    # t_stell[t_stell<0]=0
    # t_stell=t_stell[:,int(allothetic_dur):] #ignore first 100 ms
    decoded_pos=((lamb/(2*np.pi))*((np.angle(np.sum((t_stell*np.exp(1j*cell_phases[:,np.newaxis])),axis=0)))))%lamb
    decoded_pos_unwrapped=np.unwrap(decoded_pos,period=lamb)
    x= (np.arange(0,params['sim_dur'])[:])
    slope=stats.linregress(x[:params["extra_params"]["dir_change_t"]],np.unwrap(decoded_pos_unwrapped[:params["extra_params"]["dir_change_t"]])).slope
    
    #calculate positions based on speed of the network
    const_vel =slope
    init_t_r = allothetic_dur
    init_pos_r = decoded_pos[init_t_r]
    init_t_l = params["extra_params"]["dir_change_t"]+3000
    init_pos_l = decoded_pos[init_t_l]
    true_pos = (const_vel*(x-init_t_r)+init_pos_r)%(lamb)
    true_pos_l = (-const_vel*(x-init_t_l)+init_pos_l)%(lamb)
    true_pos[init_t_l:]=true_pos_l[init_t_l:]
    return true_pos

def plot_left_vs_right(stell_spks_l,params,true_pos,cell_to_analyse=92,win_size=6):
    import matplotlib.pyplot as plt
    bounds=(params["extra_params"]["dir_change_t"],params["sim_dur"])
    allothetic_dur= params["allothetic_dur"]
    cell_spks_r = np.array(stell_spks_l[cell_to_analyse])
    cell_spks_l = np.array(stell_spks_l[cell_to_analyse+params["N_per_sheet"]])
    separate_fields_r = remove_uneven_fields(cell_spks_r,(allothetic_dur,bounds[0]))
    separate_fields_l = remove_uneven_fields(cell_spks_l,(bounds[0]+allothetic_dur,bounds[1]))
    field_pos_l = calculate_field_positions(separate_fields_l,true_pos)
    field_pos_r = calculate_field_positions(separate_fields_r,true_pos)
    convoled_l = convole_field_pos(field_pos_l,params,win_size=win_size)
    convoled_l=np.mean(convoled_l,axis=0)
    convoled_l=convoled_l
    convoled_r = convole_field_pos(field_pos_r,params,win_size=win_size)
    convoled_r=np.mean(convoled_r,axis=0)
    convoled_r=convoled_r
    pos = np.arange(0,params["lambda0"],params["lambda0"]/params["n_phases"])
    means = np.mean(np.row_stack((convoled_l,convoled_r)),axis=0)
    fig,axs = plt.subplots(1,2,figsize=(12,4))
    axs[0].eventplot(field_pos_r)
    axs[0].eventplot(field_pos_l,color="red")
    axs[0].set_xlim([0,pos[-1]])


    axs[1].plot(pos,convoled_r)
    axs[1].plot(pos,convoled_l)
    axs[1].plot(pos,means)
    peak_l=pos[np.argmax(convoled_l)]
    peak_r=pos[np.argmax(convoled_r)]
    peak_mean=pos[np.argmax(means)]
    plt.show()

def circular_convolve(x, h):
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
    # return np.real(np.fft.ifft( np.fft.fft(x)*np.fft.fft(y)))

def calc_bias_stell_intrnrn(stell_decoded,intrnrn_decoded,params,bounds):
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
    l_bias=-l_bias #bias irection is opposite in other direction
    bias = np.mean(np.concatenate((r_bias,l_bias)))
    return bias
    
    
    
def calc_predictive_code(stell_spikes_l,intrnrn_spks_l,params,sim_num):
    bounds = (params["extra_params"]["dir_change_t"],int(params["sim_dur"]))
    # params["lambda0"]=64
    params["allothetic_dur"]=3000
    true_pos = decode_pos_by_intrnrn(intrnrn_spks_l,params,win_size=40)
    stell_decoded = a_utils.decode_pos(stell_spikes_l,params)
    bias=calc_bias_stell_intrnrn(stell_decoded,true_pos,params,bounds)
    print(sim_num,flush=True)
    return bias
    
def calc_inhib_g_at_first_and_last_spike(stell_spks,stell_syn_inhib_g):
    """
    stell_spks: one ring
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