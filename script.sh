#!/bin/bash

echo "Starting simulation on 18th Oct, 1:45pM"

echo "hello how are you, this is script1"
python m_sim_setup.py specs/theta_mod/NEW_THETA.py -o

echo "hello how are you, this is script2"
python m_sim_setup.py specs/theta_mod/NEW_THETA_1.py -o



echo "done"