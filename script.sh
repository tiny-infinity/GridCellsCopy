#!/bin/bash

echo "Starting simulation on 19th Oct, 12:15pM"

echo "hello how are you, this is script1"
python m_sim_setup.py specs/theta_mod/NEW_THETA_2.py -o

echo "hello how are you, this is script2"
python m_sim_setup.py specs/theta_mod/NEW_THETA_3.py -o



echo "done"