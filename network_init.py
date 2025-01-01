from neuron import h
from neuron.units import *
import numpy as np
from network import Network
import h5py
import importlib


def network_intialize(params):
    """Initialize network and setup instrumentation

    Parameters:
        params : dict
            Parameter dictionary

    Returns:
        network 
            Network object that includes cells and instrumentations

    """

    #if asked to build matrix, load from cache. Otherwise load a previously saved matrix.
    if params['build_conn_matrix']:
        file = h5py.File(
            f"cache/matrix_{params['conn_id']}_{params['sim_id']}_{params['sim_num']}.hdf5", "r"
        )
        adj_matrix = file["matrix"]
    else:
        file = h5py.File(
            f"network_configs/connections/saved_matrices/matrix_{params['conn_id']}.hdf5", "r"
        )
        adj_matrix = file["matrix"]
    network = Network(0, adj_matrix, params)  # initialize grid cells
    file.close()
    
    #Add instrumentation
    setup_instrumentation = importlib.import_module(f"network_configs.instrumentations.{params['instr_id']}_instr").setup_instrumentation
    setup_instrumentation(network)
    return network
