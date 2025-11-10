TITLE Theta Synaptic Conductance (Modulated Amplitude & Freq)

UNITS {
    (mV) = (millivolt)
    (S) = (siemens)
    (mA) = (milliamp)
    (ms) = (millisecond)
}

NEURON {
    SUFFIX thetasyn_mod
    NONSPECIFIC_CURRENT i
    RANGE g, e, i           
    GLOBAL omega, phi, dc, theta_vel_tuning, amp_vel_tuning, g_max_base
}


PARAMETER {
    g_max_base = 1e-4 (S/cm2) : The base (or max) conductance
    e = 0 (mV) : Reversal potential
    omega = 0.05 (1/ms) : Base frequency (~8 Hz)
    phi = 0
    dc = 0 : Your input variable (e.g., velocity)
    theta_vel_tuning = 0 : Flag to turn on/off freq modulation
    amp_vel_tuning = 0 : Flag to turn on/off amplitude modulation
}


ASSIGNED {
    v (mV)
    t (ms)
    i (mA/cm2)
    g (S/cm2) : The final, time-varying conductance
    freq (1/ms) : The modulated frequency
    current_g_max (S/cm2) : The modulated amplitude
}


BREAKPOINT {
    : 1. Calculate the modulated frequency
    freq = dc_to_freq(dc, theta_vel_tuning)
    
    : 2. Calculate the modulated amplitude
    current_g_max = modulate_gmax(dc, amp_vel_tuning, g_max_base)
    
    : 3. Calculate the time-varying conductance 'g'
    g = current_g_max * (1 + sin(2 * 3.14159265 * freq * (t + phi))) / 2
    
    : 4. Calculate the output current 'i'
    i = g * (v - e)
}

: --- Your frequency modulation function (unchanged) ---
FUNCTION dc_to_freq(x, y) {
    if (y == 0) {
        dc_to_freq = omega
    } else {
        dc_to_freq = (8.03) * x + 0.02556
    }
}

: --- NEW amplitude modulation function ---
FUNCTION modulate_gmax(x, y, base_gmax) {
    : x is 'dc' (your input variable)
    : y is 'amp_vel_tuning' (your flag)
    : base_gmax is 'g_max_base'
    
    if (y == 0) {
        : No modulation, just use the base amplitude
        modulate_gmax = base_gmax
    } else {
        : --- PUT YOUR MODULATION FORMULA HERE ---
        : Example: a simple linear increase with 'dc'
        : This is just a placeholder!
        : You need to replace this with your actual desired relationship.
        modulate_gmax = (0.001) * x + 0.0001 
        
        : Or, if 'base_gmax' is the maximum and you want to scale it by dc:
        : modulate_gmax = base_gmax * x 
    }
}

