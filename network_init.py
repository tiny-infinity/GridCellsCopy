from neuron import h
from neuron.units import *
import numpy as np
from network import Network
import h5py
import importlib


def network_intialize(params):
    
    if params['build_conn_matrx']:
        file = h5py.File(
            f"cache/networks/matrix_{params['connect_id']}_{params['sim_id']}_{params['sim_num']}.hdf5", "r"
        )
        adj_matrix = file["matrix"]
    else:
        file = h5py.File(
            f"network_configs/connections/saved_matrices/matrix_{params['connect_id']}.hdf5", "r"
        )
        adj_matrix = file["matrix"]
    network = Network(0, adj_matrix, params)  # initialize grid cells
    file.close()
    
    setup_instrumentation = importlib.import_module(f"network_configs/instrumentations/{params['instr_id']}_instr.py").setup_instrumentation
    setup_instrumentation(network)
    return network
