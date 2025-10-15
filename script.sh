#!/bin/bash

echo "Starting simulation on 15th Oct, 1:45AM"

echo "hello how are you, this is script1"
python m_sim_setup.py specs/theta_mod/THTA_1.py -o

echo "hello how are you, this is script2"
python m_sim_setup.py specs/theta_mod/THTA_2.py -o

echo "hello how are you, this is script3"
python m_sim_setup.py specs/theta_mod/THTA_3.py -o

echo "hello how are you, this is script4"
python m_sim_setup.py specs/theta_mod/THTA_4.py -o

echo "hello how are you, this is script5"
python m_sim_setup.py specs/theta_mod/THTA_11.py -o

echo "hello how are you, this is script6"
python m_sim_setup.py specs/theta_mod/THTA_21.py -o

echo "hello how are you, this is script5"
python m_sim_setup.py specs/theta_mod/THTA_12.py -o

echo "hello how are you, this is script6"
python m_sim_setup.py specs/theta_mod/THTA_22.py -o


echo "done"