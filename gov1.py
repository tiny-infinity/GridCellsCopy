import numpy as np
from ratinabox import Agent
from ratinabox import Environment
import matplotlib.pyplot as plt
from tqdm import tqdm

np.random.seed(42)

total_time = 30  # seconds
Env = Environment(params={
    'dimensionality':'1D',
    'boundary_conditions':'periodic',
    'scale': 2*np.pi/100   
})

Ag = Agent(Environment=Env,
           params={
               'speed_mean':0.08,   # m/s ≈ 8 cm/s
               'speed_std':0.08,
               'dt':0.00025        # s
})

for i in range(int(total_time/Ag.dt)):
    Ag.update()

pos_array = np.array(Ag.history['pos'])*100    # in cm
time_array = np.linspace(0, total_time, int(total_time/Ag.dt)) # in s

cue_pos = 3  # cm

def periodic_distance(x, y, L):
    """Shortest distance on a circle of length L"""
    return np.minimum(np.abs(x - y), L - np.abs(x - y))

def theta_amp_mod(cur_pos, cue_pos, L=2*np.pi):
    base_theta_amp = 1.0   # baseline
    min_theta_amp  = 0.1   # min at cue
    pos_var        = 10    # width of cue effect (cm)

    cur_pos = np.array(cur_pos).flatten()
    d = periodic_distance(cur_pos, cue_pos, L)

    mod_term = (base_theta_amp - min_theta_amp) * np.exp(-(d ** 2) / pos_var)
    return base_theta_amp - mod_term

def theta_oscillations(amp_array,freq,total_time,time_step):
    sine_array = np.sin(2*np.pi*freq*np.linspace(0,total_time,int(total_time/time_step)))
    return amp_array * sine_array

theta_mod_array = theta_amp_mod(pos_array, cue_pos=cue_pos)
theta_osc = theta_oscillations(theta_mod_array, freq=10, total_time=total_time, time_step=Ag.dt)

fig, axs = plt.subplots(2, 1, figsize=(8,6), sharex=False)

# trajectory in space
axs[0].plot(time_array*1000, pos_array)  # time in ms, pos in cm
axs[0].axhline(cue_pos, color='red', linestyle='--', label='Cue')
axs[0].set_xlabel("Time (ms)")
axs[0].set_ylabel("Position (cm)")
axs[0].legend()

# theta modulation in time
axs[1].plot(time_array*1000, theta_osc, color='blue')
axs[1].set_xlabel("Time (ms)")
axs[1].set_ylabel("Theta amplitude")

plt.tight_layout()
plt.show()
