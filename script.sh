#!/bin/bash

echo "Starting simulation on 13th Nov, 12:46AM"

echo "hello how are you, this is script1"
python m_sim_setup.py specs/RA_2/1e3.py -o

echo "hello how are you, this is script1"
python m_sim_setup.py specs/RA_1/7e3.py -o

echo "hello how are you, this is script1"
python m_sim_setup.py specs/RA_2/7e3.py -o

echo "hello how are you, this is script1"
python m_sim_setup.py specs/RA_1/6e3.py -o


echo "done"