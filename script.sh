#!/bin/bash

echo "Starting simulation on 19th Oct, 12:15pM"

echo "hello how are you, this is script1"
python m_sim_setup.py specs/theta_mod_s2/THETA_CONST.py -o

echo "hello how are you, this is script1"
python m_sim_setup.py specs/theta_mod_s2/THETA_MOD.py -o

echo "hello how are you, this is script1"
python m_sim_setup.py specs/theta_mod_s2/THETA_OSC.py -o




echo "done"