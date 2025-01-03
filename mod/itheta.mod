TITLE Theta for interneurons



UNITS {
    (mV)=(millivolt)
    (S) = (siemens)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX i_theta
    NONSPECIFIC_CURRENT itheta 
    RANGE Amp,dc,c
    GLOBAL vthresh,omega,phi,theta_vel_tuning
    

      
}

PARAMETER {
    Amp = 1e-4 (S/cm2)
    vthresh = -80 (mV)
    omega = 0.01 (1/ms)
    phi=0
    dc=0
    c=0
    theta_vel_tuning=0}   

ASSIGNED {
        v (mV)
        itheta (mA/cm2)

} 




BREAKPOINT { 

	itheta = Amp*sin(2*3.14*dc_to_freq(dc,theta_vel_tuning)*(t+phi)) + c
    

}

FUNCTION dc_to_freq(x,y) { 
  if (y==0) {
    dc_to_freq=omega
  } else {
        dc_to_freq=(8.03)*x+0.02556
}
}