TITLE Theta for stellate


UNITS {
    (mV)=(millivolt)
    (S) = (siemens)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX i_theta_stell
    NONSPECIFIC_CURRENT itheta 
    RANGE Amp,dc,c,omega
    GLOBAL phi
    

      
}

PARAMETER {
    Amp = 0 (S/cm2)
    vthresh = -80 (mV)
    omega = 0.01 (1/ms)
    phi=0
    c=0}   

ASSIGNED {
        v (mV)
        itheta (mA/cm2)

} 




BREAKPOINT { 

  itheta = (Amp*sin(2*3.14*omega*(t+phi)) + c)
    

}